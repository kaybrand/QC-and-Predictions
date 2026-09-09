# Local patches to installed packages

These patch files live in the repo because the things they patch do **not** —
they are edits to packages inside a conda env
(`$FORGE/envs/run_snakemake9`). `pip install -U`, `conda update`, or rebuilding
the env reverts them **silently**, and the failure they prevent does not
reappear until days into a long run. So: after any change to that env, re-check
this directory.

## `snakemake-executor-plugin-slurm-2.3.1-nodefail-exclude.patch`

**Applies to:** `snakemake-executor-plugin-slurm` **2.3.1** (pypi), in
`$FORGE/envs/run_snakemake9/lib/python3.11/site-packages/snakemake_executor_plugin_slurm/__init__.py`

**Applied:** 2026-09-09. Pristine copy kept beside the file as
`__init__.py.orig-2.3.1-preNodeFailPatch` (md5 `dab23d0cea61e482f68a552bd68f2772`).

**Revert:** `cp __init__.py.orig-2.3.1-preNodeFailPatch __init__.py`

### What it fixes

One `NODE_FAIL` can kill an entire multi-day run, hours after the fact.

On `NODE_FAIL` the plugin records the bad node so later jobs avoid it:

```python
node = sacct_output.strip()          # sacct -j <id> -n -X -o nodelist%-256
if node:
    self._failed_nodes.add(node)
```

When a job never got a node, `sacct` returns not a hostname but the sentinel
**`None assigned`**. It is non-empty, so it is cached as though it were a node.
Every later submission then gets, unquoted:

```
sbatch ... --exclude=None assigned --wrap="..."
```

`sbatch` parses `assigned` as a positional script filename and dies with
`sbatch: error: Unable to open file assigned`. Nothing can be submitted for the
rest of the run.

**This happened:** run `42023901` (catlas-human, 258 clusters) ran cleanly from
2026-09-04 11:00 to 2026-09-07 00:00, then three `NODE_FAIL`s poisoned the set
and it spent ~6 h failing every submission before giving up at 2d19h of a 7-day
allocation. 2045 of its 2747 errors were that one message, and the only
`--exclude=` value in the whole 16 MB log was `None assigned`. It produced zero
new final prediction sets.

Note the value is *transient*: `sacct` now reports real node names for all three
of those jobs. Do not expect to reproduce it by re-querying after the fact.

### The two changes

1. **Quote the value** (`shlex.quote`, already imported). Empirically the
   sufficient fix — verified with `sbatch --test-only`: `--exclude='None
   assigned'` is **accepted** and schedules normally, `--exclude=None assigned`
   reproduces the error exactly. Slurm tolerates an unknown node name in
   `--exclude`; the outage was purely the shell word-split. Also protects the
   bracket-range nodelist form (`sh04-16n[25-27]`), whose `[...]` is a shell
   glob class.

2. **Reject the sentinel** before it is cached. Rejects on *whitespace* rather
   than the literal string, because a real Slurm NodeList never contains any —
   so the whole sentinel class is covered, not one example. Without this the
   junk entry survives for the life of the run, every submission carries a
   pointless `--exclude`, and the shutdown warning misreports `None assigned` as
   a failed node, in exactly the log you would read during the next incident.

Verified after applying: module imports, `snakemake --version` still 9.16.3, and
the guard rejects `"None assigned" / "None" / "none" / "" / "   "` while
admitting `sh04-16n25`, `sh04-16n[25-27]`, `sh04-16n25,sh04-16n26`.

### Upstream

**Still unfixed in 2.8.0, the latest release** (checked 2026-09-08 at tag
`v2.8.0`). The code moved — the `--exclude` construction now lives in
`submit_string.py:157`, the admission in `__init__.py:1526` — but both are
byte-for-byte the same logic. So **upgrading does not fix this**, and there is
no reason to upgrade for it.

Upgrading also would not require a newer Snakemake, if it ever becomes
worthwhile: the plugin declares no dependency on `snakemake` at all, only on the
interface packages, and even 2.8.0's floors (`snakemake-interface-common>=1.23.0`,
`snakemake-interface-executor-plugins>=9.4.0`) sit inside the ranges snakemake
9.16.3 itself allows (`>=1.20.1,<2.0` and `>=9.3.2,<10.0`). Versions 2.4.0–2.5.4
need no dependency changes at all; 2.6.x needs `jobstep>=0.6.0`.

Worth a PR. `submit_string.py` already imports a `safe_quote` helper and uses it
on the neighbouring `--reservation`, `--qos`, `-C` and `-D` — `--exclude` is the
one that was missed, so the fix has an in-file precedent.
