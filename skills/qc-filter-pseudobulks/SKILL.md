---
name: qc-filter-pseudobulks
description: Make or redo per-cluster QC guides from human-chosen thresholds for IGVF multiome pseudobulks, and filter them to QC-passing ATAC fragments and RNA count matrices. Routes to the qc-and-predictions skill, which now covers this as Steps 2–5 of the full pipeline.
allowed-tools: Read, Write, Edit, Bash, Glob, Grep
---

# QC Filter Pseudobulks (moved)

This workflow is now part of **`qc-and-predictions`**, which runs the QC-and-Predictions
pipeline end to end on a new IGVF primary-pseudobulk dataset.

1. Read `../qc-and-predictions/SKILL.md` first. It holds the hard rules (never upload, never
   expose credentials, writable paths) that apply to this work too.
2. For QC guides, thresholds, merged clusters, the pipeline config and Cell Annotation
   overrides, follow `../qc-and-predictions/references/qc-guides.md`.

Kept under its old name so existing `/qc-filter-pseudobulks` invocations still land
somewhere useful.
