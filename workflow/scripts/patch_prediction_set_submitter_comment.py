#!/usr/bin/env python
"""One-off backfill: adds a "Version 1" marker to submitter_comment on every
already-live PredictionSet on the portal (e.g. IGVF4's, which predate
prediction_set.py's own constant_fields["submitter_comment"] = "Version 1" --
see that table module -- and so never got the marker at POST time).

Deliberately NOT folded into orchestrator.plan_table's normal POST/PATCH
flow: that flow always recomputes a row's full payload from local config
alone (orchestrator.build_payload), with no awareness of a field's CURRENT
live value -- fine for every other field, wrong here, since blindly
re-patching submitter_comment to the constant "Version 1" on every
reconciliation run would clobber any genuine free text a submitter_comment
already carries (out-of-band portal edits, or a future per-cluster comment
analogous to filtered_barcode_list.py's cluster_cfg.get("submitter_comment")).

Idempotent tag-and-check rule, applied to each PredictionSet's CURRENT
submitter_comment (fetched live, never assumed from local state):
  - empty/None            -> "Version 1"
  - already contains "Version 1" (plain substring) -> unchanged, skipped
  - anything else         -> "(Version 1) - {existing text}"
Re-running this script (this migration again, or a future one reusing the
same check) always lands on the second branch and no-ops -- the marker never
stacks ("(Version 1) - (Version 1) - ..." can't happen) and a submitter's
own free text is only ever prefixed once, never overwritten.

Matches on the literal substring "Version 1", not a stricter anchored tag
-- accepted tradeoff (see PredictionSet submitter_comment design discussion):
a submitter_comment that already happens to mention "Version 1" for an
unrelated reason (e.g. "compare against Version 1 of the other pipeline")
would be mistaken for an already-applied marker and left alone. Considered
acceptable since submitter_comment is free text nobody currently populates
this way, and a stricter anchor (e.g. requiring the literal prefix
"(Version 1) - ") would then fail to recognize genuinely already-marked
comments that predate this convention.
"""

import argparse
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from igvf_metadata import portal_client  # noqa: E402

VERSION_MARKER = "Version 1"
PROFILE_ID = "prediction_set"

_QUERY = (
    "type=PredictionSet&status%21=deleted&limit=all"
    "&field=%40id&field=aliases&field=submitter_comment&field=uuid&field=accession"
)


def log(msg):
    print(f"[patch_prediction_set_submitter_comment] {msg}", file=sys.stderr)


def _alias_of(row):
    aliases = row.get("aliases") or []
    return aliases[0] if aliases else None


def _dataset_of(alias):
    """alias is "{alias_prefix}:{dataset}_{cluster}_scE2G_{family}_predictions"
    (context.make_alias joins parts with "_") -- dataset is always a bare
    igvfN token with no internal underscore, so the first "_"-delimited
    piece after the ":" reliably separates it out."""
    if alias is None or ":" not in alias:
        return None
    suffix = alias.split(":", 1)[1]
    return suffix.split("_", 1)[0]


def _record_id(row):
    return row.get("uuid") or row.get("accession") or row.get("@id")


def new_submitter_comment(current):
    """Returns the value to PATCH to, or None if no change is needed."""
    if not current:
        return VERSION_MARKER
    if VERSION_MARKER in current:
        return None
    return f"({VERSION_MARKER}) - {current}"


def plan(rows, datasets):
    to_patch = []
    counts = {"skipped-no-alias": 0, "skipped-wrong-dataset": 0, "unchanged": 0, "set-fresh": 0, "prepended": 0}
    for row in rows:
        alias = _alias_of(row)
        if alias is None:
            counts["skipped-no-alias"] += 1
            log(f"WARNING: PredictionSet {row.get('@id')} has no alias -- skipping")
            continue
        dataset = _dataset_of(alias)
        if datasets and dataset not in datasets:
            counts["skipped-wrong-dataset"] += 1
            continue
        current = row.get("submitter_comment")
        new = new_submitter_comment(current)
        if new is None:
            counts["unchanged"] += 1
            continue
        counts["set-fresh" if not current else "prepended"] += 1
        to_patch.append(
            {"record_id": _record_id(row), "alias": alias, "old": current, "new": new}
        )
    return to_patch, counts


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument(
        "--dataset", default="",
        help='comma-separated dataset token(s) to restrict to (e.g. "igvf4"); default: every live PredictionSet',
    )
    p.add_argument("-o", "--output", default="prediction_set_submitter_comment_patch.tsv", help="patch TSV path")
    p.add_argument(
        "--mode", choices=["preview", "validate", "upload"], default="preview",
        help="preview (default): write the patch TSV only. validate: also run iu_register.py --dry-run. "
        "upload: PATCH for real -- requires explicitly choosing this.",
    )
    p.add_argument("--igvf-mode", default="prod", help="passed through to igvf_utils -- defaults to prod")
    p.add_argument("--iu-register-path", default=portal_client.IU_REGISTER_DEFAULT_PATH)
    args = p.parse_args()

    datasets = {d for d in args.dataset.split(",") if d}
    reader = portal_client.PortalReader(igvf_mode=args.igvf_mode)

    log("fetching live PredictionSet objects...")
    rows = reader.get_multireport(_QUERY)
    log(f"{len(rows)} PredictionSet(s) fetched" + (f", filtering to dataset(s) {sorted(datasets)}" if datasets else ""))

    to_patch, counts = plan(rows, datasets)
    log(f"plan: {counts}, {len(to_patch)} to patch")
    for entry in to_patch:
        log(f"  {entry['alias']}: {entry['old']!r} -> {entry['new']!r}")

    written = portal_client.write_tsv(
        args.output,
        [{"submitter_comment": e["new"]} for e in to_patch],
        record_ids=[e["record_id"] for e in to_patch],
    )
    if written is None:
        log("nothing to patch -- no file written")
        return
    log(f"wrote {written}")

    if args.mode == "preview":
        return

    result = portal_client.invoke_register(
        written, PROFILE_ID, patch=True, dry_run=(args.mode == "validate"),
        igvf_mode=args.igvf_mode, iu_register_path=args.iu_register_path,
    )
    log(f"iu_register.py exit={result.returncode}")
    if result.returncode != 0:
        log(result.stderr[-2000:])
        sys.exit(result.returncode)

    if args.mode == "upload":
        mismatches = 0
        for entry in to_patch:
            # database=True: a same-process re-GET against the default
            # Elasticsearch-backed index can lag a just-completed PATCH by
            # several seconds and falsely report the pre-patch value --
            # confirmed 2026-08-11 against real production data. Read the
            # database directly instead.
            record = reader.get_by_alias(entry["alias"], database=True)
            live_comment = record.get("submitter_comment") if record else None
            if live_comment != entry["new"]:
                mismatches += 1
                log(f"MISMATCH after patch: {entry['alias']}: expected {entry['new']!r}, portal has {live_comment!r}")
        log(f"verified {len(to_patch) - mismatches}/{len(to_patch)} patched submitter_comment(s)")


if __name__ == "__main__":
    main()
