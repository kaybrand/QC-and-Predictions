"""Shared config/context objects for the IGVF metadata upload backbone.

Every table module builds its rows from one Context per (dataset, cluster,
model-or-None) scope-key, so table code never has to know how IgvfConfig
defaults were resolved or where scE2G's outputs live on disk.
"""

import os
from dataclasses import dataclass, field
from typing import Optional

# This repo's own root -- mirrors workflow/rules/common.smk's
# `WDIR = os.path.dirname(workflow.basedir)`, computed the equivalent way
# since this package has no access to Snakemake's `workflow` object.
WDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))

# The `main`-branch QC_pseudobulks worktree's own root -- a sibling worktree
# of this repo (this repo is the igvf-portal-submission worktree), per `git
# worktree list`. multiome_data_cluster_dir below reads already-filtered
# ATAC/RNA data from there, same "read pre-existing artifacts from the
# main-branch worktree" pattern as tables/qc_documents.py's
# QC_PSEUDOBULKS_PLOTS_DIR (fixed 2026-08-03 after the same class of bug:
# this used to be `WDIR`-relative, silently finding nothing since this
# worktree has no multiome_data/ dir of its own).
QC_PSEUDOBULKS_WDIR = "/oak/stanford/groups/engreitz/Projects/IGVF-E2GPillarProject/QC_pseudobulks"


@dataclass(frozen=True)
class IgvfConfig:
    """lab/award/alias_prefix are identical across all ~20 tables. Defaults
    are this lab's real values; override via the `igvf:` block in a
    *_pipeline_config.yaml so a different lab can reuse these scripts
    without editing code.

    enabled_families gates which scE2G model families' cluster_model-scoped
    rows (Prediction Tabular Files, Signal Files, BEDPE Index File,
    Prediction Set) actually get generated this run -- 2026-07-20 feedback:
    only Multiome is uploaded this year even though a cluster's own `models`
    config may list scATAC too (that list reflects what scE2G ran, not what
    IGVF should receive). Enforced centrally in orchestrator._iter_scopes."""

    lab: str = "/labs/jesse-engreitz/"
    award: str = "/awards/HG011972/"
    alias_prefix: str = "jesse-engreitz"
    enabled_families: tuple = ("Multiome",)

    @classmethod
    def from_dict(cls, d):
        d = d or {}
        return cls(
            lab=d.get("lab", cls.lab),
            award=d.get("award", cls.award),
            alias_prefix=d.get("alias_prefix", cls.alias_prefix),
            enabled_families=tuple(d.get("enabled_families", cls.enabled_families)),
        )


@dataclass
class Context:
    dataset: str
    cluster: str
    model: Optional[str]  # None for cluster-scoped tables (scope="cluster")
    cluster_cfg: dict  # config["clusters"][dataset][cluster]: models, pseudobulk_annotation, qc_guide
    igvf: IgvfConfig
    scE2G_dir: str
    cache: dict = field(default_factory=dict)  # per-run memoization (e.g. score thresholds, one lookup per model)
    conn: Optional[object] = None  # state.db connection -- cell_metadata.get_metadata_for's cache lookup needs it

    @property
    def results_dir(self):
        return os.path.join(self.scE2G_dir, "results", "uniformly_processed")

    @property
    def cluster_dir(self):
        return os.path.join(self.results_dir, self.dataset, self.cluster)

    @property
    def multiome_data_cluster_dir(self):
        """The Synapse-side filtered_data location, NOT scE2G's own results
        dir. Corrected 2026-08-03: reads from QC_PSEUDOBULKS_WDIR (the
        main-branch QC_pseudobulks worktree), not this worktree's own WDIR --
        this worktree has no multiome_data/ dir of its own; the real
        already-filtered ATAC/RNA data for existing clusters lives in the
        main-branch worktree (produced by that repo's legacy manual filtering
        before this Snakemake pipeline existed).

        NOTE this now DIVERGES from workflow/rules/common.smk's own
        multiome_data_dir(dataset), which is still WDIR-relative (this
        worktree) -- that's what filter_pseudobulks.smk's rules actually
        write fresh output to. So a cluster with no pre-existing legacy data
        that gets filtered for the first time by THIS pipeline's own Snakemake
        rules would land in a different directory than this property looks
        in. Flagged, not resolved -- unifying the two (or teaching this
        property to check both locations) is unresolved follow-up."""
        return os.path.join(QC_PSEUDOBULKS_WDIR, "multiome_data", self.dataset, self.cluster)

    def with_model(self, model):
        return Context(
            self.dataset, self.cluster, model, self.cluster_cfg, self.igvf, self.scE2G_dir, self.cache, self.conn
        )


def make_alias(igvf: IgvfConfig, *parts) -> str:
    """Every ITEM_ALIAS across every table is "{alias_prefix}:{'_'.join(parts)}" --
    this is the one place that format is defined.

    UNRESOLVED, flagged 2026-08-03 for follow-up in the coming months (do
    NOT change `dataset` pre-emptively -- nothing here is confirmed yet):
    this pipeline's informal `dataset` labels (igvf1, igvf2, ...) are
    expected to be retired from IGVF Portal-facing identifiers sometime in
    2027, in favor of each dataset's principal analysis set accession (e.g.
    igvf4 -> IGVFDS5875AFXS). Every alias built here -- and the Kundaje-lab
    "anshul-kundaje:{dataset}-{cluster}-{subsample}"-shaped aliases in
    refs.py/tables/principal_pseudobulk_set.py -- embeds `dataset` literally,
    so all of them are exposed once that rename lands. The current
    dataset<->accession mapping is regenerable via
    igvf_cell_annotation_report/map_dataset_to_principal_analysis_set/
    check_dataset_accession_mapping.py (writes
    dataset_to_principal_analysis_set_accession.json in that same
    directory) -- useful for scoping the eventual migration, not needed for
    anything today."""
    return f"{igvf.alias_prefix}:" + "_".join(str(p) for p in parts)
