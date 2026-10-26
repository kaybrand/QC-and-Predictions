# QC-and-Predictions — Claude skills branch

This orphan branch holds only Claude Code skills for this repo, versioned
independently of the code branches. It is checked out once, in its own worktree
(`IGVF-E2GPillarProject/QC-and-Predictions-skills/`), and never switched or
merged into a code branch.

- `skills/qc-and-predictions/` — run the pipeline end to end on a new IGVF
  primary-pseudobulk dataset; `references/qc-guides.md` has the QC-guide detail.
- `skills/qc-filter-pseudobulks/` — pointer kept for the old skill name.

Claude loads these through the EngreitzLabAgents plugin: the local checkout's
`skills/<name>` is a git-ignored symlink to `skills/<name>` here, so every session
(any worktree, any branch, any directory) uses this one copy. Edit and commit here.
