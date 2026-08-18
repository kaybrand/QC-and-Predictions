"""Shared classification of IGVF Portal PseudobulkSet records.

Moved here from cell_metadata.py (2026-08-17) so the two consumers -- that
module's Cell Annotation cache and portal_files.py's download discovery --
cannot drift apart on the one rule that decides what a pseudobulk set IS.
Kept deliberately free of any dependency on the rest of this package, so
importing it can never introduce a cycle.

The rule, confirmed against real production data:

  primary   -- input_file_sets are all AnalysisSets. Produced upstream (the
               Kundaje lab) directly from analysis sets. These are the
               pipeline's INPUTS.
  principal -- input_file_sets are all PseudobulkSets. Produced by THIS
               pipeline, derived from primaries. These are our own OUTPUTS.
  None      -- anything else: no input_file_sets at all, or a mix.

Two traps this encodes, both verified 2026-08-17 against api.data.igvf.org:

1. **"Primary" is not a label on the portal.** Unlike analysis sets (which
   carry file_set_type "principal analysis" / "intermediate analysis"), a
   primary pseudobulk set is indistinguishable from a principal one by
   file_set_type alone -- BOTH are "pseudobulk analysis" (all 1069 primary and
   all 14 principal sets share that value). The input_file_sets @id prefix is
   the only discriminator. Do not "simplify" this to a file_set_type check.

2. **input_file_sets entries carry no "@type" sub-field**, even when
   requested, so classification must go by the "@id" PATH PREFIX. An earlier
   draft of the caller assumed "@type" was available and silently classified
   everything as None.

Real distribution of the 1710 non-deleted PseudobulkSets (2026-08-17):
  primary 1069 | mixed 627 | principal 14
The 627 "mixed" are all-/curated-sets/ CATlas mouse pseudobulks
(yang-li:catlas-mouse-pseudobulk-*, input file_set_type "external sequencing
data") -- correctly excluded here; they are a separate workstream with a
different filename convention. A dataset appearing in that bucket is expected,
not a bug, but callers should COUNT and report it rather than dropping it
silently.
"""

ANALYSIS_SET_PREFIX = "/analysis-sets/"
PSEUDOBULK_SET_PREFIX = "/pseudobulk-sets/"


def input_file_set_ids(row):
    """The "@id" of every input_file_sets entry, skipping anything that isn't a
    dict (the portal embeds these as objects; a bare string would mean the
    caller forgot field=input_file_sets)."""
    return [e.get("@id", "") for e in row.get("input_file_sets") or [] if isinstance(e, dict)]


def classify(row):
    """"primary", "principal", or None -- see this module's docstring."""
    ids = input_file_set_ids(row)
    if not ids:
        return None
    if all(i.startswith(ANALYSIS_SET_PREFIX) for i in ids):
        return "primary"
    if all(i.startswith(PSEUDOBULK_SET_PREFIX) for i in ids):
        return "principal"
    return None


def input_file_set_collections(row):
    """The distinct portal collection of each input, e.g. {"analysis-sets"} or
    {"curated-sets"} -- for REPORTING which bucket an unclassifiable row fell
    into, so a growing population of skipped sets is visible instead of silent."""
    return sorted({i.strip("/").split("/")[0] for i in input_file_set_ids(row) if i})


def principal_analysis_set(row):
    """The accession of this set's `principal analysis` input, if it has
    exactly one -- e.g. IGVFDS5875AFXS for an igvf4 pseudobulk, matching
    igvf_cell_annotation_report/map_dataset_to_principal_analysis_set/
    dataset_to_principal_analysis_set_accession.json.

    This is the identity `dataset` is expected to become (see
    context.make_alias's UNRESOLVED note), so it is worth recording now even
    though nothing reads it yet. Returns None rather than guessing when a set
    has zero or several -- a real primary set also carries `intermediate
    analysis` inputs (the example IGVFDS0317IYCP has one principal and two
    intermediate), so this must filter on file_set_type, not just take the
    first input."""
    accessions = {
        e.get("accession")
        for e in row.get("input_file_sets") or []
        if isinstance(e, dict) and e.get("file_set_type") == "principal analysis" and e.get("accession")
    }
    return next(iter(accessions)) if len(accessions) == 1 else None
