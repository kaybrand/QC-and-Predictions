"""Discovery of the pseudobulk component files to download from the IGVF Portal.

One read-only multireport GET (via portal_client.PortalReader.get_multireport,
unchanged) enumerates every PseudobulkSet with its member `files` embedded, then
this module selects the PRIMARY sets' ATAC/RNA/per-cell-QC files and works out
where each one belongs on local disk.

Everything below was verified against api.data.igvf.org on 2026-08-17 with an
authenticated request. Counts are from that run (1710 non-deleted PseudobulkSets).

WHAT WE SELECT, and why each filter exists
------------------------------------------
1. PRIMARY sets only (pseudobulk_sets.classify == "primary": input_file_sets are
   all AnalysisSets). 1069 of 1710. PRINCIPAL sets (14) are *this pipeline's own
   uploads* derived from primaries -- downloading them would re-fetch our own
   output. The remaining 627 are CATlas mouse pseudobulks under /curated-sets/,
   a separate workstream. See pseudobulk_sets.py for why file_set_type cannot
   make this distinction.

2. One lab (default DEFAULT_LAB, Anshul Kundaje's). Of the 1069 primary sets,
   1031 are Kundaje's and hold ALL 3087 target files; the other 38 -- from four
   other labs, including chongyuan-luo:pgp_start_C29_* -- carry NONE of the
   three target content_types. So this filter and "has target files" coincide
   exactly today; it is here to keep a new lab's differently-shaped pseudobulks
   from silently entering the pipeline's input tree.

3. content_type in TARGET_CONTENT_TYPES. Deliberately NOT also filtered on
   file_format: across all 1031 sets each target content_type has exactly one
   format and one filename, with zero variation, so a format filter would add no
   selectivity -- and if upstream ever DOES change a format, that is a reason to
   notice (EXPECTED_FILENAMES mismatch -> needs_review) rather than to skip the
   file and silently lose a pseudobulk.

4. NOT filtered on `status` being released. Every target file is "in progress",
   which is a normal internal state (it is what necessitates the "Version 1"
   submitter comments), not a quality signal -- filtering to `released` would
   download nothing at all. Only genuinely superseded states are dropped
   (EXCLUDED_FILE_STATUSES); `deleted` is already excluded by the query itself.

WHERE EACH FILE GOES
--------------------
`submitted_file_name` is the path oracle: 3084 of 3087 target files match
_ABS_SFN below, from a single root, and it is the only field carrying dataset,
annotation, subsample AND filename together. Cross-checked against the set alias
("anshul-kundaje:{dataset}-{annotation}-{subsample}"), the two agreed on all 3084
-- so the alias is used purely as a corroborating check, never as the source.

The pseudobulk DIRECTORY NAME is reused VERBATIM from that path rather than
rebuilt from parsed parts. Rebuilding would normalise away exactly the oddities
worth preserving (see the annotation_0 case below) and risks a subtle mismatch
with the archive we compare against, where directory identity is the join key.

Deliberately NOT used: `samples`. 34 primary sets have more than one sample, so
samples[0] is not a safe subsample source -- and submitted_file_name makes it
unnecessary. Also deliberately not used: `href` and Content-Disposition, for the
filename -- the S3 object is named "{accession}.{ext}", not the pipeline filename.

KNOWN IRREGULAR SETS (do not "fix" by skipping)
-----------------------------------------------
IGVFDS1981DGID / IGVFDS3848ZSPZ / IGVFDS7151WSNC, aliased
"anshul-kundaje:pseudobulk-IGVFDS8428QVAO-endothelial cell-IGVFSM...", are real
primary pseudobulks in an in-progress state: each currently has ONLY a
`cell by gene matrix` (plus gene quantifications) and is missing `fragments` and
`per-cell quality report`, and its submitted_file_name is RELATIVE
("pseudobulks/annotation_0-IGVFSM0220TLGT/rna_counts_mtx.h5ad") with no dataset
component. They are expected to gain the missing files. So: download what exists,
resolve the dataset from the set's principal analysis set accession, flag
needs_review, and let a later run complete them automatically. Note their aliases
also embed a SPACE and the accession-style dataset identity that `dataset` is
slated to become -- another reason never to parse identity out of an alias here.
"""

import os
import re
import sys
import urllib.parse
from collections import Counter

from . import pseudobulk_sets

DEFAULT_LAB = "/labs/anshul-kundaje/"

# content_type -> the filename upstream uses for it. The filename is taken from
# submitted_file_name, never from this map; the map is the EXPECTATION, used only
# to flag a file whose name has changed shape (-> needs_review).
TARGET_CONTENT_TYPES = {
    "fragments": "fragments.tsv.gz",
    "cell by gene matrix": "rna_counts_mtx.h5ad",
    "per-cell quality report": "per_cell_qc.tsv.gz",
}

# Superseded/withdrawn files. "deleted" is already excluded by the query's
# status!=deleted. "in progress" is NOT here on purpose -- see module docstring.
EXCLUDED_FILE_STATUSES = frozenset({"revoked", "replaced"})

_MULTIREPORT_FIELDS = (
    "@id",
    "accession",
    "aliases",
    "status",
    "file_set_type",
    "lab",
    "input_file_sets",
    "samples",
    "files",
    "cell_annotation",
    "cell_qualifier",
    "cell_type",
)

# Absolute form, 3084/3087 of target files, single root
# (/oak/stanford/groups/kasowski/sbaek1/igvf/). `dirname` is captured whole and
# reused verbatim -- see module docstring.
_ABS_SFN = re.compile(
    r"^.*/igvf/(?P<dataset>[^/]+)/pseudobulks/(?P<dirname>[^/]+)/(?P<filename>[^/]+)$"
)
# Relative fallback, the 3 in-progress sets above: no dataset component.
_REL_SFN = re.compile(r"^(?:.*/)?pseudobulks/(?P<dirname>[^/]+)/(?P<filename>[^/]+)$")
# A well-formed pseudobulk directory. Anchored on the IGVFSM subsample accession
# at the end, so a hyphen inside the annotation name is harmless -- never a
# blind split on "-".
_DIRNAME = re.compile(r"^annotation-(?P<annotation>.+)-(?P<subsample>IGVFSM\w+)$")


def log(msg):
    print(f"[portal_files] {msg}", file=sys.stderr)


def multireport_query():
    """type=PseudobulkSet with files embedded. `field=files` alone returns fully
    embedded file objects (same as cell_metadata's input_file_sets/samples/
    cell_type) -- no sub-field enumeration needed. limit=all is required and
    verified sufficient: it returned all 1710 rows in one request, and no
    pagination mechanism exists anywhere in igvf_utils to fall back on."""
    params = [("type", "PseudobulkSet"), ("status!", "deleted"), ("limit", "all")]
    params += [("field", f) for f in _MULTIREPORT_FIELDS]
    return urllib.parse.urlencode(params)


def _alias_of(row):
    aliases = row.get("aliases") or []
    return aliases[0] if aliases else None


def resolve_scope(row, file_obj):
    """Where this file belongs on disk, as a path RELATIVE to the download root:
    "{dataset}/pseudobulks/{dirname}/{filename}".

    Returns (scope_dict, review_reasons). scope_dict always carries
    submitted_file_name; dataset/annotation/subsample/dirname/filename/rel_path
    are None when they could not be determined WITHOUT guessing. review_reasons
    is a list of short tags -- non-empty means a human should look, and the
    caller marks the row needs_review rather than dropping it."""
    sfn = file_obj.get("submitted_file_name") or ""
    reasons = []
    scope = {
        "submitted_file_name": sfn or None,
        "dataset": None,
        "annotation": None,
        "subsample": None,
        "dirname": None,
        "filename": None,
        "rel_path": None,
    }
    if not sfn:
        reasons.append("no_submitted_file_name")
        return scope, reasons

    m = _ABS_SFN.match(sfn)
    if m:
        scope["dataset"] = m.group("dataset")
    else:
        m = _REL_SFN.match(sfn)
        if not m:
            reasons.append("unparseable_submitted_file_name")
            return scope, reasons
        reasons.append("relative_submitted_file_name")
        # No dataset in the path. The principal analysis set accession is the
        # only non-guessed identity available, and is what `dataset` is slated
        # to become anyway.
        pas = pseudobulk_sets.principal_analysis_set(row)
        if pas:
            scope["dataset"] = pas
            reasons.append("dataset_from_principal_analysis_set")
        else:
            reasons.append("no_dataset")

    scope["dirname"] = m.group("dirname")
    scope["filename"] = m.group("filename")

    dm = _DIRNAME.match(scope["dirname"])
    if dm:
        scope["annotation"] = dm.group("annotation")
        scope["subsample"] = dm.group("subsample")
    else:
        # e.g. "annotation_0-IGVFSM0220TLGT" -- keep the directory verbatim and
        # say so, rather than inventing an annotation name.
        reasons.append("nonstandard_directory_name")

    expected = TARGET_CONTENT_TYPES.get(file_obj.get("content_type"))
    if expected and scope["filename"] != expected:
        # Not fatal: a changed upstream format is something to surface, not to
        # drop. The real filename is still what gets written.
        reasons.append(f"unexpected_filename:{scope['filename']}")

    alias = _alias_of(row)
    if alias and scope["subsample"] and scope["subsample"] not in alias:
        reasons.append("alias_subsample_mismatch")

    if scope["dataset"] and scope["dirname"] and scope["filename"]:
        scope["rel_path"] = os.path.join(
            scope["dataset"], "pseudobulks", scope["dirname"], scope["filename"]
        )
    return scope, reasons


def discover(reader, lab=DEFAULT_LAB, content_types=None, datasets=None):
    """Returns (records, report).

    records: one dict per selected file, ready for state.upsert_portal_file plus
    a "rel_path" and "review_reasons" the caller resolves against its own
    download root. Nothing is dropped for being irregular -- irregular rows come
    back with review_reasons set and rel_path possibly None.

    report: counters describing everything seen AND everything excluded, so no
    exclusion is silent.
    """
    wanted = set(content_types or TARGET_CONTENT_TYPES)
    unknown = wanted - set(TARGET_CONTENT_TYPES)
    if unknown:
        raise ValueError(f"unknown content_type(s): {sorted(unknown)}")

    rows = reader.get_multireport(multireport_query())
    log(f"multireport returned {len(rows)} PseudobulkSet row(s)")

    report = {
        "sets_total": len(rows),
        "sets_by_class": Counter(),
        "unclassified_by_collection": Counter(),
        "primary_by_lab": Counter(),
        "files_by_content_type": Counter(),
        "files_excluded_by_status": Counter(),
        "sets_missing_content_type": Counter(),
        "review_reasons": Counter(),
        "datasets": Counter(),
        "sets_selected": 0,
        "sets_with_no_target_files": 0,
    }
    records = []

    for row in rows:
        kind = pseudobulk_sets.classify(row)
        report["sets_by_class"][kind or "unclassified"] += 1
        if kind is None:
            report["unclassified_by_collection"][
                ",".join(pseudobulk_sets.input_file_set_collections(row)) or "none"
            ] += 1
            continue
        if kind != "primary":
            continue

        row_lab = (row.get("lab") or {}).get("@id")
        report["primary_by_lab"][row_lab] += 1
        if lab and row_lab != lab:
            continue

        pas = pseudobulk_sets.principal_analysis_set(row)
        alias = _alias_of(row)
        set_accession = row.get("accession")
        seen_types = set()
        set_records = []

        for f in row.get("files") or []:
            ct = f.get("content_type")
            if ct not in wanted:
                continue
            status = f.get("status")
            if status in EXCLUDED_FILE_STATUSES:
                report["files_excluded_by_status"][f"{ct}:{status}"] += 1
                continue
            seen_types.add(ct)
            scope, reasons = resolve_scope(row, f)
            for r in reasons:
                report["review_reasons"][r.split(":", 1)[0]] += 1
            report["files_by_content_type"][ct] += 1
            if scope["dataset"]:
                report["datasets"][scope["dataset"]] += 1
            set_records.append(
                {
                    "accession": f.get("accession"),
                    "file_set": row.get("@id"),
                    "principal_analysis_set": pas,
                    "lab": row_lab,
                    "content_type": ct,
                    "file_format": f.get("file_format"),
                    "href": f.get("href"),
                    "alias": (f.get("aliases") or [None])[0],
                    "md5sum": f.get("md5sum"),
                    "file_size": f.get("file_size"),
                    "portal_status": status,
                    "upload_status": f.get("upload_status"),
                    "submitted_file_name": scope["submitted_file_name"],
                    "dataset": scope["dataset"],
                    "annotation": scope["annotation"],
                    "subsample": scope["subsample"],
                    "rel_path": scope["rel_path"],
                    "review_reasons": reasons,
                    "set_accession": set_accession,
                    "set_alias": alias,
                }
            )

        # A set restricted away by --datasets is not "missing" anything; only
        # judge completeness for sets we actually kept.
        if datasets is not None:
            set_records = [r for r in set_records if r["dataset"] in datasets]
            if not set_records:
                continue

        if not set_records:
            report["sets_with_no_target_files"] += 1
            continue

        missing = wanted - seen_types
        if missing:
            # Expected and temporary for the in-progress sets described in the
            # module docstring: report it, keep what exists, let a later run
            # complete the set.
            report["sets_missing_content_type"][",".join(sorted(missing))] += 1
            log(
                f"INCOMPLETE {set_accession} ({alias}): has {sorted(seen_types)}, "
                f"missing {sorted(missing)} -- keeping what exists"
            )
            for r in set_records:
                r["review_reasons"] = r["review_reasons"] + [
                    "set_missing:" + ",".join(sorted(missing))
                ]

        report["sets_selected"] += 1
        records.extend(set_records)

    return records, report


def log_report(report):
    log(f"sets by classification: {dict(report['sets_by_class'])}")
    if report["unclassified_by_collection"]:
        log(f"  unclassified sets by input collection: {dict(report['unclassified_by_collection'])}")
    log(f"primary sets by lab: {dict(report['primary_by_lab'])}")
    log(f"primary sets selected: {report['sets_selected']}")
    log(f"files by content_type: {dict(report['files_by_content_type'])}")
    if report["files_excluded_by_status"]:
        log(f"  files excluded by status: {dict(report['files_excluded_by_status'])}")
    if report["sets_with_no_target_files"]:
        log(f"  sets carrying none of the target content_types: {report['sets_with_no_target_files']}")
    if report["sets_missing_content_type"]:
        log(f"  sets missing some target content_type: {dict(report['sets_missing_content_type'])}")
    if report["review_reasons"]:
        log(f"  needs-review reasons: {dict(report['review_reasons'])}")
    log(f"files per dataset: {dict(sorted(report['datasets'].items()))}")
