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


@dataclass(frozen=True)
class IgvfConfig:
    """lab/award/alias_prefix are identical across all ~20 tables. Defaults
    are this lab's real values; override via the `igvf:` block in a
    *_pipeline_config.yaml so a different lab can reuse these scripts
    without editing code."""

    lab: str = "/labs/jesse-engreitz/"
    award: str = "/awards/HG011972"
    alias_prefix: str = "jesse-engreitz"

    @classmethod
    def from_dict(cls, d):
        d = d or {}
        return cls(
            lab=d.get("lab", cls.lab),
            award=d.get("award", cls.award),
            alias_prefix=d.get("alias_prefix", cls.alias_prefix),
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

    @property
    def results_dir(self):
        return os.path.join(self.scE2G_dir, "results", "uniformly_processed")

    @property
    def cluster_dir(self):
        return os.path.join(self.results_dir, self.dataset, self.cluster)

    @property
    def multiome_data_cluster_dir(self):
        """Mirrors common.smk's multiome_data_dir(dataset)/cluster -- the
        Synapse-side filtered_data location, NOT scE2G's own results dir."""
        return os.path.join(WDIR, "multiome_data", self.dataset, self.cluster)

    def with_model(self, model):
        return Context(self.dataset, self.cluster, model, self.cluster_cfg, self.igvf, self.scE2G_dir, self.cache)


def make_alias(igvf: IgvfConfig, *parts) -> str:
    """Every ITEM_ALIAS across every table is "{alias_prefix}:{'_'.join(parts)}" --
    this is the one place that format is defined."""
    return f"{igvf.alias_prefix}:" + "_".join(str(p) for p in parts)
