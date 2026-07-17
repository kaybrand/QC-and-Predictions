"""Declarative registry for the ~20 IGVF metadata tables.

Each table is one TableSpec: a portal object type, a scope (one row-group
per cluster, or per cluster+model), a handful of field layers merged in a
fixed order, and a list of VariantSpecs -- the actual row-producing units
within that table (e.g. "full"/"thresholded"/"bedpe"/"elements"/"genes" for
Prediction Tabular Files).

Field layering (see orchestrator.build_payload):
  base fields (aliases/award/lab -- identical for every row, every table)
  -> table.constant_fields (identical for every row of THIS table)
  -> table.scope_fields(ctx) (identical across variants for one scope-key,
     e.g. file_set/cell_type_annotation)
  -> variant.build_row(ctx) (the part that actually varies row to row)

A row's dependencies (variant.depends_on) name other (table_name, variant)
pairs that must already be status='uploaded' in the state ledger. Because
run() just gets re-invoked (once per pipeline-triggered cluster, or
periodically by the scanner), a row whose dependency isn't ready yet is
simply left 'pending' -- no explicit topological sort needed, only that
depends_on correctly names what has to exist first.
"""

from dataclasses import dataclass, field
from typing import Callable, List, Tuple

Dependency = Tuple[str, str]  # (table_name, variant_name) that must be status='uploaded' first

_REGISTRY = {}


@dataclass
class VariantSpec:
    name: str
    build_row: Callable  # (ctx) -> dict[str, Any]
    enabled: Callable = lambda ctx: True  # (ctx) -> bool -- should this row exist at all this run
    depends_on: Callable = lambda ctx: []  # (ctx) -> list[Dependency]


@dataclass
class TableSpec:
    name: str  # registry key, e.g. "prediction_tabular_files"
    object_type: str  # portal profile/schema name, e.g. "tabular_file"
    scope: str  # "cluster" | "cluster_model"
    build_alias: Callable  # (ctx, variant_name) -> str
    required_columns: List[str]  # columns WE must populate (schema-required but blank-ok, e.g. md5sum, excluded)
    variants: List[VariantSpec]
    constant_fields: dict = field(default_factory=dict)
    scope_fields: Callable = lambda ctx: {}  # (ctx) -> dict, same across all variants of one scope-key


def register(spec: TableSpec) -> TableSpec:
    if spec.name in _REGISTRY:
        raise ValueError(f"table {spec.name!r} already registered")
    _REGISTRY[spec.name] = spec
    return spec


def get(name: str) -> TableSpec:
    return _REGISTRY[name]


def all_specs() -> List[TableSpec]:
    return list(_REGISTRY.values())
