"""Discovery of the pseudobulk component files to download from the IGVF Portal.

One read-only multireport GET (via portal_client.PortalReader.get_multireport,
unchanged) enumerates every PseudobulkSet with its member `files` embedded, then
this module works out which files to fetch and where each one belongs on disk.

TWO FAMILIES, TWO PATH ORACLES
------------------------------
`discover()`        -- the PRIMARY sets' ATAC/RNA/per-cell-QC files (Kundaje
                       lab). Path oracle is `submitted_file_name`; the set
                       alias is only a corroborating cross-check.
`discover_catlas()` -- the CATlas curated-set pseudobulks (Yang Li lab), one
                       `fragments` file each. Path oracle is the SET ALIAS;
                       `submitted_file_name` is the corroborating cross-check.
                       That inversion is deliberate and load-bearing -- see
                       that function's own docstring.

The two never share a run: they select disjoint sets (`classify` returns
"primary" for one and None for the other), and each has its own default scope,
so neither can silently pull in the other's pseudobulks.

Everything in the `discover()` section below was verified against
api.data.igvf.org on 2026-08-17 with an authenticated request; counts are from
that run (1710 non-deleted PseudobulkSets, since grown to 1849 as of
2026-09-02 -- see pseudobulk_sets.py, and do not treat any total here as an
invariant). The CATlas section was verified live on 2026-09-02.

WHAT WE SELECT, and why each filter exists
------------------------------------------
1. PRIMARY sets only (pseudobulk_sets.classify == "primary": input_file_sets are
   all AnalysisSets). 1069 of 1710. PRINCIPAL sets (14) are *this pipeline's own
   uploads* derived from primaries -- downloading them would re-fetch our own
   output. The remaining 627 are CATlas pseudobulks under /curated-sets/ --
   369 mouse AND 258 human, NOT mouse-only as this line used to claim -- handled
   by discover_catlas() below, not by this function. See pseudobulk_sets.py for
   why file_set_type cannot make this distinction.

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
`per-cell quality report`. They are expected to gain the missing files, so:
download what exists, flag needs_review, and let a later run complete them.

Three separate irregularities in those three sets, each handled explicitly:

  - submitted_file_name is RELATIVE with no dataset component
    ("pseudobulks/annotation_0-IGVFSM0220TLGT/rna_counts_mtx.h5ad"), so `dataset`
    comes from the set's principal analysis set accession (IGVFDS8428QVAO).
  - "annotation_0" IS NOT AN ANNOTATION. It is the submitter's positional
    placeholder (annotation index 0) and appears verbatim in the portal's own
    submitted_file_name. The real annotation is "endothelial cell", taken from
    the SET ALIAS via annotation_from_alias() and independently corroborated by
    cell_type.term_name (CL:0000115). An earlier version of this module recorded
    "annotation_0-IGVFSM..." as the annotation, which was simply wrong.
  - The alias uses the future "{accession}-{term name}-{subsample}" shape and the
    annotation contains a SPACE. That is why annotation_from_alias peels from both
    ends using the already-known dataset and subsample instead of splitting on
    "-": a positional split cannot survive a value with spaces or hyphens.

Note the on-disk path deliberately keeps the VERBATIM upstream directory name
("annotation_0-IGVFSM0220TLGT"), not the recovered annotation. Two reasons: the
mirror stays faithful to where upstream actually put the bytes, and it keeps a
space out of a filesystem path. The recovered annotation is carried in the
`annotation` field, which is what the comparison and manifest group on -- so
"what this pseudobulk IS" and "where its bytes live" stay separate.

Also worth knowing before comparing them to anything: these three are Mus
musculus, Parse SPLiT-seq ("Parse SPLiT-seq pseudobulk of Mus musculus adrenal
gland endothelial cell"), unlike every other set here, which is human 10x
multiome.
"""

import os
import re
import sys
import urllib.parse
from collections import Counter

from . import pseudobulk_sets

DEFAULT_LAB = "/labs/anshul-kundaje/"

# --- CATlas (Yang Li lab) curated-set pseudobulks -------------------------
CATLAS_LAB = "/labs/yang-li/"

# "yang-li:catlas-{species}-pseudobulk-{cluster}".
#
# `cluster` is (?P<cluster>.+) -- deliberately greedy-to-end and NEVER split on
# "-", because real cluster names contain hyphens: "class-06_CTX-CGE_GABA",
# "156_MB-ant-ve_Dmrta2_Glut", "supercluster-Bergmann_glia". `species` is
# [a-z]+, which cannot match a hyphen, so the FIRST "-pseudobulk-" after
# "catlas-{lowercase}" anchors the split unambiguously.
#
# Verified live 2026-09-02 against all 627 CATlas sets: every set has exactly
# ONE alias, all 627 match this pattern, `species` is only "human"/"mouse", and
# (species, cluster) is unique -- 0 collisions. No non-conforming set exists
# today, which is precisely why a non-match must stay LOUD (needs_review) rather
# than being silently normalised.
CATLAS_ALIAS = re.compile(r"^yang-li:catlas-(?P<species>[a-z]+)-pseudobulk-(?P<cluster>.+)$")

# The species this pipeline recognises. A new species appearing upstream must be
# an explicit decision (it needs its own genome build, annotations and manifest
# reference_files), never an automatic inclusion -- so an unknown species token
# is reported, not processed.
CATLAS_SPECIES = ("human", "mouse")

# CATlas sets also carry "annotated cell by peak matrix" (hdf5), "cell
# annotations" (tsv) and "signal of unique reads" (bigWig). scE2G consumes only
# the fragments, so only fragments is a target; the others are deliberately not
# fetched. The expected filename is a FORMAT STRING, not a fixed name, because
# CATlas embeds species+cluster in it -- see catlas_expected_filename.
CATLAS_TARGET_CONTENT_TYPES = ("fragments",)

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
# Looser: ANY directory ending in a subsample accession. Recovers the subsample
# from a directory whose annotation part is not a real annotation name -- e.g.
# "annotation_0-IGVFSM0220TLGT", where "annotation_0" is the submitter's
# positional placeholder (annotation index 0), not a cell type.
_DIRNAME_SUBSAMPLE_ONLY = re.compile(r"^.*?-(?P<subsample>IGVFSM\w+)$")


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


def annotation_from_alias(alias, dataset, subsample):
    """Recover the annotation (cluster) from a SET alias when the directory name
    doesn't carry a usable one.

    The alias is "{lab}:[pseudobulk-]{dataset}-{annotation}-{subsample}", e.g.
    "anshul-kundaje:pseudobulk-IGVFDS8428QVAO-endothelial cell-IGVFSM0220TLGT"
    -> annotation "endothelial cell".

    Peeled from BOTH ENDS using values already known independently (the dataset
    and the subsample accession), never by splitting on "-" and taking a
    position. That matters here: the annotation itself contains a space, and in
    the general case may contain hyphens, so any positional split is wrong.

    Returns None rather than a guess if the alias doesn't have the expected
    shape. Cross-check only -- callers must not let this override an annotation
    the directory name already gave them."""
    if not alias or not subsample:
        return None
    suffix = alias.split(":", 1)[-1]
    tail = f"-{subsample}"
    if not suffix.endswith(tail):
        return None
    suffix = suffix[: -len(tail)]
    # Optional literal segment seen on the accession-style aliases.
    if suffix.startswith("pseudobulk-"):
        suffix = suffix[len("pseudobulk-") :]
    if dataset and suffix.startswith(f"{dataset}-"):
        suffix = suffix[len(dataset) + 1 :]
    return suffix or None


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

    alias = _alias_of(row)
    dm = _DIRNAME.match(scope["dirname"])
    if dm:
        scope["annotation"] = dm.group("annotation")
        scope["subsample"] = dm.group("subsample")
    else:
        # The directory's annotation segment is not an annotation. Real case:
        # "annotation_0-IGVFSM0220TLGT", where "annotation_0" is the submitter's
        # positional placeholder (annotation index 0). The subsample accession at
        # the end is still trustworthy, and the ANNOTATION is recoverable from the
        # set alias -- for these sets it is "endothelial cell", corroborated by
        # cell_type.term_name (CL:0000115).
        reasons.append("nonstandard_directory_name")
        sm = _DIRNAME_SUBSAMPLE_ONLY.match(scope["dirname"])
        if sm:
            scope["subsample"] = sm.group("subsample")
        else:
            samples = [x.get("accession") for x in (row.get("samples") or []) if isinstance(x, dict)]
            if len(samples) == 1:
                scope["subsample"] = samples[0]
        annotation = annotation_from_alias(alias, scope["dataset"], scope["subsample"])
        if annotation:
            scope["annotation"] = annotation
            reasons.append("annotation_from_alias")
        else:
            reasons.append("no_annotation")

    expected = TARGET_CONTENT_TYPES.get(file_obj.get("content_type"))
    if expected and scope["filename"] != expected:
        # Not fatal: a changed upstream format is something to surface, not to
        # drop. The real filename is still what gets written.
        reasons.append(f"unexpected_filename:{scope['filename']}")

    if alias and scope["subsample"] and scope["subsample"] not in alias:
        reasons.append("alias_subsample_mismatch")

    if scope["dataset"] and scope["dirname"] and scope["filename"]:
        scope["rel_path"] = os.path.join(
            scope["dataset"], "pseudobulks", scope["dirname"], scope["filename"]
        )
    return scope, reasons


def parse_catlas_alias(alias):
    """(species, cluster) from a CATlas set alias, or (None, None) if it doesn't
    conform. See CATLAS_ALIAS for why this must not split on "-"."""
    if not alias:
        return None, None
    m = CATLAS_ALIAS.match(alias)
    if not m:
        return None, None
    return m.group("species"), m.group("cluster")


def catlas_dataset(species):
    """The local dataset label for a CATlas species: "catlas-human".

    This is the seam that makes species mixing structurally impossible rather
    than merely discouraged: every path, config block and Portal alias is keyed
    on this label, so one dataset can only ever contain one species."""
    return f"catlas-{species}"


def catlas_expected_filename(species, cluster):
    """What upstream names a CATlas fragments file. Verified to match all 627
    sets' submitted_file_name basename exactly (2026-09-02). Used ONLY as a
    cross-check against the alias, never as the source of species/cluster."""
    return f"{species}-{cluster}_fragments.tsv.gz"


def catlas_resolve_scope(row, file_obj, species, cluster):
    """Where a CATlas fragments file belongs, relative to the download root:
    "{catlas-species}/fragments/{filename}".

    THE PATH ORACLE IS THE SET ALIAS, not submitted_file_name -- the inverse of
    resolve_scope() above, and deliberately so. CATlas's submitted_file_name is
    a relative, non-pseudobulk-layout path ("0.file/human-Bgl_51_fragments.tsv.gz",
    verified on all 627): it has no `pseudobulks/` component and no dataset
    component, so neither _ABS_SFN nor _REL_SFN matches it and it cannot yield a
    dataset at all. The alias can, and does so for 627/627 with zero exceptions.

    submitted_file_name is still used for two things: the FILENAME (taken
    verbatim, same discipline as the Kundaje path -- mirror what upstream
    actually named the bytes) and as a CROSS-CHECK, since it independently
    encodes species+cluster. A disagreement between the two means one of our two
    assumptions broke, which is a needs_review, never a silent pick-one.

    There is no annotation/subsample dir level here: a CATlas set is one cluster
    with one fragments file, so the Kundaje layout's
    "annotation-{annotation}-{subsample}" directory has nothing to represent.
    `subsample` stays None; `annotation` carries the cluster, which is what the
    manifest and comparison layers group on."""
    sfn = file_obj.get("submitted_file_name") or ""
    reasons = []
    scope = {
        "submitted_file_name": sfn or None,
        "dataset": None,
        "annotation": cluster,
        "subsample": None,
        "dirname": None,
        "filename": None,
        "rel_path": None,
    }
    if not species or not cluster:
        reasons.append("unparseable_catlas_alias")
        return scope, reasons
    if species not in CATLAS_SPECIES:
        # Not fatal, but must not be processed as if we knew what to do with it:
        # a new species needs its own genome build and manifest reference_files.
        reasons.append(f"unknown_species:{species}")

    scope["dataset"] = catlas_dataset(species)
    expected = catlas_expected_filename(species, cluster)
    if sfn:
        scope["filename"] = os.path.basename(sfn)
        if scope["filename"] != expected:
            reasons.append(f"unexpected_filename:{scope['filename']}")
    else:
        reasons.append("no_submitted_file_name")
        scope["filename"] = expected
        reasons.append("filename_from_alias")

    upload_status = file_obj.get("upload_status")
    if upload_status != "validated":
        # 3 of the 627 are "invalidated" as of 2026-09-02 (the portal's own
        # checkfiles rejected them), one of which is human
        # (supercluster-Bergmann_glia). Surfaced per-row so the ledger is
        # self-describing, but still ATTEMPTED -- same "download what exists,
        # flag it" stance the irregular Kundaje sets get. The md5 check is what
        # actually decides whether the bytes are usable.
        reasons.append(f"upload_status:{upload_status}")

    scope["rel_path"] = os.path.join(scope["dataset"], "fragments", scope["filename"])
    return scope, reasons


def _file_record(row, file_obj, scope, reasons, row_lab, set_accession, set_alias, pas):
    """The one definition of a discovery record's shape. Both discover() and
    discover_catlas() build rows through here so the two families cannot drift
    apart on the keys download_portal_pseudobulks.record_discovery and
    state.upsert_portal_file require."""
    return {
        "accession": file_obj.get("accession"),
        "file_set": row.get("@id"),
        "principal_analysis_set": pas,
        "lab": row_lab,
        "content_type": file_obj.get("content_type"),
        "file_format": file_obj.get("file_format"),
        "href": file_obj.get("href"),
        "s3_uri": file_obj.get("s3_uri"),
        "alias": (file_obj.get("aliases") or [None])[0],
        "md5sum": file_obj.get("md5sum"),
        "file_size": file_obj.get("file_size"),
        "portal_status": file_obj.get("status"),
        "upload_status": file_obj.get("upload_status"),
        "submitted_file_name": scope["submitted_file_name"],
        "dataset": scope["dataset"],
        "annotation": scope["annotation"],
        "subsample": scope["subsample"],
        "rel_path": scope["rel_path"],
        "review_reasons": reasons,
        "set_accession": set_accession,
        "set_alias": set_alias,
    }


def _empty_report(total):
    """The report shape both discover() and discover_catlas() fill, so
    log_report() can print either. Counters describe everything seen AND
    everything excluded, so no exclusion is silent."""
    return {
        "sets_total": total,
        "sets_by_class": Counter(),
        "unclassified_by_collection": Counter(),
        "primary_by_lab": Counter(),
        "files_by_content_type": Counter(),
        "files_excluded_by_status": Counter(),
        # Surfaced because it is the only field that says anything about whether
        # the BYTES exist -- href and s3_uri are calculated from record identity
        # and are present even for a file whose S3 object 404s (verified across
        # all 9832 files on every PseudobulkSet: both fields present 100% of the
        # time, for every upload_status value). "validated" is the healthy value;
        # "invalidated" means the portal's own checkfiles rejected the file.
        # 477 Kundaje `peaks` files are in that state, and 3 of the 627 CATlas
        # fragments files -- out of scope to fix here, but a TARGET content_type
        # landing there must never pass silently.
        "upload_status": Counter(),
        "sets_missing_content_type": Counter(),
        "review_reasons": Counter(),
        "datasets": Counter(),
        "sets_selected": 0,
        "sets_with_no_target_files": 0,
    }


def discover_catlas(reader, species=None, content_types=None):
    """Returns (records, report) for the CATlas (Yang Li lab) curated-set
    pseudobulks. Same contract as discover(): nothing is dropped for being
    irregular, irregular rows come back with review_reasons set.

    species: None (THE DEFAULT) selects EVERY CATlas cluster, both species --
        627 sets as of 2026-09-02. Passing "human" (or an iterable of species)
        selects ONLY those, and is the supported way to run one species.

    Species can never mix within a dataset: the dataset label is
    "catlas-{species}" (catlas_dataset), so the species is part of the identity
    of every path and alias downstream rather than a filter applied over a
    shared namespace. Restricting by species narrows which datasets appear; it
    can never produce a dataset holding two species.

    A CATlas set is skipped ONLY for being a different species than requested,
    or for not being a CATlas set at all. Every other irregularity is reported.
    """
    wanted_species = None
    if species is not None:
        wanted_species = {species} if isinstance(species, str) else set(species)
        unknown = wanted_species - set(CATLAS_SPECIES)
        if unknown:
            raise ValueError(
                f"unknown CATlas species {sorted(unknown)}; known: {list(CATLAS_SPECIES)}"
            )

    wanted_types = set(content_types or CATLAS_TARGET_CONTENT_TYPES)
    unknown_types = wanted_types - set(CATLAS_TARGET_CONTENT_TYPES)
    if unknown_types:
        raise ValueError(
            f"unknown CATlas content_type(s): {sorted(unknown_types)}; "
            f"known: {list(CATLAS_TARGET_CONTENT_TYPES)}"
        )

    rows = reader.get_multireport(multireport_query())
    log(f"multireport returned {len(rows)} PseudobulkSet row(s)")

    report = _empty_report(len(rows))
    report["catlas_sets_by_species"] = Counter()
    report["catlas_sets_skipped_by_species"] = Counter()
    report["catlas_nonconforming_aliases"] = []
    records = []

    for row in rows:
        row_lab = (row.get("lab") or {}).get("@id")
        if row_lab != CATLAS_LAB:
            continue
        alias = _alias_of(row)
        set_species, cluster = parse_catlas_alias(alias)
        if set_species is None:
            # Zero of these exist today. If one ever appears it is a genuine
            # upstream convention change and must be loud, because the alias is
            # this family's only path oracle -- there is no fallback that could
            # place the file correctly.
            report["catlas_nonconforming_aliases"].append(
                {"set": row.get("accession"), "alias": alias}
            )
            log(
                f"NON-CONFORMING CATlas alias on {row.get('accession')}: {alias!r} -- "
                "cannot resolve species/cluster, skipping (see portal_files.CATLAS_ALIAS)"
            )
            continue

        report["catlas_sets_by_species"][set_species] += 1
        if wanted_species is not None and set_species not in wanted_species:
            report["catlas_sets_skipped_by_species"][set_species] += 1
            continue

        report["sets_by_class"][pseudobulk_sets.classify(row) or "unclassified"] += 1
        pas = pseudobulk_sets.principal_analysis_set(row)
        set_accession = row.get("accession")
        seen_types, set_records, excluded_by_status = set(), [], []

        for f in row.get("files") or []:
            ct = f.get("content_type")
            if ct not in wanted_types:
                continue
            status = f.get("status")
            if status in EXCLUDED_FILE_STATUSES:
                excluded_by_status.append(f"{ct}:{status}")
                continue
            seen_types.add(ct)
            scope, reasons = catlas_resolve_scope(row, f, set_species, cluster)
            set_records.append(
                _file_record(row, f, scope, reasons, row_lab, set_accession, alias, pas)
            )

        if not set_records:
            report["sets_with_no_target_files"] += 1
            continue

        missing = wanted_types - seen_types
        if missing:
            report["sets_missing_content_type"][",".join(sorted(missing))] += 1
            log(
                f"INCOMPLETE {set_accession} ({alias}): has {sorted(seen_types)}, "
                f"missing {sorted(missing)} -- keeping what exists"
            )
            for r in set_records:
                r["review_reasons"] = r["review_reasons"] + [
                    "set_missing:" + ",".join(sorted(missing))
                ]

        for rec in set_records:
            report["files_by_content_type"][rec["content_type"]] += 1
            report["upload_status"][rec["upload_status"]] += 1
            if rec["dataset"]:
                report["datasets"][rec["dataset"]] += 1
            for reason in rec["review_reasons"]:
                report["review_reasons"][reason.split(":", 1)[0]] += 1
        for tag in excluded_by_status:
            report["files_excluded_by_status"][tag] += 1

        report["sets_selected"] += 1
        records.extend(set_records)

    return records, report


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

    report = _empty_report(len(rows))
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

        excluded_by_status = []
        for f in row.get("files") or []:
            ct = f.get("content_type")
            if ct not in wanted:
                continue
            status = f.get("status")
            if status in EXCLUDED_FILE_STATUSES:
                excluded_by_status.append(f"{ct}:{status}")
                continue
            seen_types.add(ct)
            scope, reasons = resolve_scope(row, f)
            set_records.append(
                _file_record(row, f, scope, reasons, row_lab, set_accession, alias, pas)
            )

        # The `datasets` restriction is applied BEFORE any per-file accounting
        # below, so every counter in `report` describes what this invocation
        # actually selected. Counting first would make a --datasets igvf10 run
        # report 13 sets alongside all 3087 files, which reads as a bug.
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

        # Accounted for LAST, so every counter reflects exactly the records this
        # invocation returns -- after the --datasets restriction, and after the
        # set_missing tag above has been attached.
        for rec in set_records:
            report["files_by_content_type"][rec["content_type"]] += 1
            report["upload_status"][rec["upload_status"]] += 1
            if rec["dataset"]:
                report["datasets"][rec["dataset"]] += 1
            for reason in rec["review_reasons"]:
                report["review_reasons"][reason.split(":", 1)[0]] += 1
        for tag in excluded_by_status:
            report["files_excluded_by_status"][tag] += 1

        report["sets_selected"] += 1
        records.extend(set_records)

    return records, report


def log_report(report):
    log(f"sets by classification: {dict(report['sets_by_class'])}")
    if report["unclassified_by_collection"]:
        log(f"  unclassified sets by input collection: {dict(report['unclassified_by_collection'])}")
    if report["primary_by_lab"]:
        log(f"primary sets by lab: {dict(report['primary_by_lab'])}")
    # CATlas-only keys, absent from a discover() report -- .get() so one
    # function can print either family without branching on which ran.
    by_species = report.get("catlas_sets_by_species")
    if by_species:
        log(f"CATlas sets by species: {dict(by_species)}")
    skipped_species = report.get("catlas_sets_skipped_by_species")
    if skipped_species:
        log(
            f"  skipped as a non-requested species: {dict(skipped_species)} "
            "(this is the --species filter working, not a problem)"
        )
    nonconforming = report.get("catlas_nonconforming_aliases")
    if nonconforming:
        log(
            f"  WARNING: {len(nonconforming)} CATlas set(s) have a NON-CONFORMING alias and were "
            "skipped -- the alias is this family's only path oracle, so these cannot be placed:"
        )
        for item in nonconforming[:10]:
            log(f"    {item['set']}: {item['alias']!r}")
    log(f"sets selected: {report['sets_selected']}")
    log(f"files by content_type: {dict(report['files_by_content_type'])}")
    upload = dict(report["upload_status"])
    log(f"files by portal upload_status: {upload}")
    unhealthy = {k: v for k, v in upload.items() if k != "validated"}
    if unhealthy:
        log(
            f"  WARNING: {sum(unhealthy.values())} selected file(s) are not upload_status "
            f"'validated': {unhealthy} -- 'file not found' will 404, 'invalidated' means the "
            "portal's checkfiles rejected the upload"
        )
    if report["files_excluded_by_status"]:
        log(f"  files excluded by status: {dict(report['files_excluded_by_status'])}")
    if report["sets_with_no_target_files"]:
        log(f"  sets carrying none of the target content_types: {report['sets_with_no_target_files']}")
    if report["sets_missing_content_type"]:
        log(f"  sets missing some target content_type: {dict(report['sets_missing_content_type'])}")
    if report["review_reasons"]:
        log(f"  needs-review reasons: {dict(report['review_reasons'])}")
    log(f"files per dataset: {dict(sorted(report['datasets'].items()))}")
