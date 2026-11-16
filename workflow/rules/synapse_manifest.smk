"""
Synapse manifest generation/upload for all four shared result spaces,
across every dataset touched in this run. Each rule shells out to
manage_synapse_manifest.py, which does the actual live Synapse querying,
scoping, diffing, and (conditionally) syncing/deleting -- see that script's
docstring for the full mechanism. Cluster identity throughout is the
(dataset, cluster) pair, passed as "dataset/cluster" tokens.

Everything here defaults to dry-run (`synapse.dry_run: true`) and to NOT
deleting/overwriting (`synapse.confirm_delete`/`confirm_overwrite: false`),
matching this pipeline's "review before pushing to shared infrastructure"
posture. These rules genuinely hit the network at DAG-execution time (unlike
the parse-time config writers in common.smk), so they only run when actually
targeted, not on every `snakemake -n`.

Each manifest output is wrapped in Snakemake's built-in `update()` flag.
Without it, Snakemake unconditionally deletes a rule's declared output file
right before running its shell command (`Job.remove_existing_output` --
confirmed by reading snakemake/jobs.py directly), which would erase this
file's pre-existing content before manage_synapse_manifest.py ever gets a
chance to read and merge it -- silently defeating that script's own
merge-with-existing-content logic every single time, regardless of how
correct that logic is. `update()` is Snakemake's own supported mechanism for
"this output gets modified in place, don't wipe it first."
"""

SYNAPSE_CFG = config.get("synapse", {})
DRY_RUN = str(SYNAPSE_CFG.get("dry_run", True))
CONFIRM_DELETE = str(SYNAPSE_CFG.get("confirm_delete", False))
CONFIRM_OVERWRITE = str(SYNAPSE_CFG.get("confirm_overwrite", False))
# Opt-in (2026-11-16): `rule all` only includes Synapse manifest targets when
# synapse.enabled is true, so a pipeline run never needs a Synapse connection
# (manage_synapse_manifest.py logs in even under --dry-run true). The rules stay
# runnable by name either way.
SYNAPSE_ENABLED = bool(SYNAPSE_CFG.get("enabled", False))
MANIFEST_DIR = os.path.join(RESULTS_DIR_BASE, "synapse_manifests")
MODEL_TOKENS = "multiome_powerlaw_v3,scATAC_powerlaw_v3"
# Reuses --model-tokens' substring-match mechanism to scope
# manage_element_list_manifest's live "owned" query to ONLY element-list
# entities in the shared predictions_parent_id folder -- "element_list"
# appears in every {dataset}_{cluster}_element_list.bed(.gz)(.tbi) filename
# and in no e2g prediction/gene-list filename, so this can't accidentally
# match into the sibling predictions manifest's own files, and vice versa.
ELEMENT_LIST_TOKENS = "element_list"


def cluster_keys_arg():
    return ",".join(f"{dataset}/{cluster}" for dataset, cluster in UPLOAD_ELIGIBLE_CLUSTERS)


def get_filtered_data_files():
    files = []
    for dataset, cluster in UPLOAD_ELIGIBLE_CLUSTERS:
        cluster_dir = os.path.join(multiome_data_dir(dataset), cluster)
        files += [
            os.path.join(cluster_dir, f"atac_fragments_{dataset}_{cluster}.tsv.gz"),
            os.path.join(cluster_dir, f"atac_fragments_{dataset}_{cluster}.tsv.gz.tbi"),
        ]
        if config["clusters"][dataset][cluster]["models"] != ["scATAC_powerlaw_v3"]:
            rna_dir = os.path.join(cluster_dir, f"rna_count_matrix_{dataset}_{cluster}")
            files += [
                os.path.join(rna_dir, "matrix.mtx.gz"),
                os.path.join(rna_dir, "barcodes.tsv.gz"),
                os.path.join(rna_dir, "features.tsv.gz"),
            ]
    return files


def get_manifest_targets():
    return [
        os.path.join(MANIFEST_DIR, "filtered_data_manifest.tsv"),
        os.path.join(MANIFEST_DIR, "predictions_manifest.tsv"),
        os.path.join(MANIFEST_DIR, "candidates_manifest.tsv"),
        os.path.join(MANIFEST_DIR, "features_manifest.tsv"),
    ]


def get_predictions_manifest_target():
    """Re-hooks manage_predictions_manifest into `rule all` (dropped from the
    default target 2026-08-03 -- see the Snakefile's own docstring). Scoped to
    just the predictions manifest, not all four of get_manifest_targets(): that
    is the one covering get_reformat_prediction_files() (the reformatted,
    portal-format scE2G prediction sets, thresholded and unthresholded -- NOT
    element lists, see get_element_list_manifest_target() below), which is
    what needs to stay in sync with Synapse as clusters finish reformatting.

    Empty unless synapse.enabled is true (see SYNAPSE_ENABLED), and empty in
    local_only mode (matching get_reformat_prediction_files() itself
    returning [] there) and whenever there is nothing upload-eligible yet --
    local_only's whole point is zero network contact, and manage_synapse_manifest.py
    logs in to Synapse unconditionally, even under --dry-run true, the moment it
    runs at all."""
    if not SYNAPSE_ENABLED or PIPELINE_MODE == "local_only" or not UPLOAD_ELIGIBLE_CLUSTERS:
        return []
    return [os.path.join(MANIFEST_DIR, "predictions_manifest.tsv")]


def get_element_list_manifest_target():
    """Same gating as get_predictions_manifest_target(), for
    manage_element_list_manifest -- split into its own manifest 2026-09-21 so
    element lists can be reviewed/synced independently of prediction sets
    (they land in the same Synapse parent folder, just a different local
    review artifact and a separate --model-tokens scope, see
    ELEMENT_LIST_TOKENS above)."""
    if not SYNAPSE_ENABLED or PIPELINE_MODE == "local_only" or not UPLOAD_ELIGIBLE_CLUSTERS:
        return []
    return [os.path.join(MANIFEST_DIR, "element_list_manifest.tsv")]


rule manage_filtered_data_manifest:
    input:
        get_filtered_data_files(),
    output:
        update(os.path.join(MANIFEST_DIR, "filtered_data_manifest.tsv")),
    params:
        parent_id=SYNAPSE_CFG.get("filtered_data_parent_id", ""),
        cluster_keys=cluster_keys_arg(),
    conda:
        "../envs/post_process.yaml"
    resources:
        mem_mb=determine_mem_mb,
    shell:
        """
        python {workflow.basedir}/scripts/manage_synapse_manifest.py \
        --product filtered_data \
        --parent-id {params.parent_id} \
        --cluster-keys "{params.cluster_keys}" \
        --nested \
        --manifest-out {output} \
        --dry-run {DRY_RUN} \
        --confirm-delete {CONFIRM_DELETE} \
        --confirm-overwrite {CONFIRM_OVERWRITE} \
        {input}
        """


rule manage_predictions_manifest:
    input:
        get_reformat_prediction_files(),
    output:
        update(os.path.join(MANIFEST_DIR, "predictions_manifest.tsv")),
    params:
        parent_id=SYNAPSE_CFG.get("predictions_parent_id", ""),
        cluster_keys=cluster_keys_arg(),
    conda:
        "../envs/post_process.yaml"
    resources:
        mem_mb=determine_mem_mb,
    shell:
        """
        python {workflow.basedir}/scripts/manage_synapse_manifest.py \
        --product predictions \
        --parent-id {params.parent_id} \
        --cluster-keys "{params.cluster_keys}" \
        --nested \
        --model-tokens {MODEL_TOKENS} \
        --manifest-out {output} \
        --dry-run {DRY_RUN} \
        --confirm-delete {CONFIRM_DELETE} \
        --confirm-overwrite {CONFIRM_OVERWRITE} \
        {input}
        """


rule manage_element_list_manifest:
    # Same Synapse parent folder as manage_predictions_manifest (element lists
    # are nested under the same {dataset}/{cluster}/ as the e2g files) -- kept
    # as its own rule/manifest/review-artifact so element lists sync on their
    # own schedule, independent of prediction-set changes. ELEMENT_LIST_TOKENS
    # scopes the live "owned" query to just element-list entities, so this
    # rule can never see (and therefore never flag as stale) the sibling
    # predictions manifest's own files in the same folder, or vice versa.
    input:
        get_reformat_element_list_files(),
    output:
        update(os.path.join(MANIFEST_DIR, "element_list_manifest.tsv")),
    params:
        parent_id=SYNAPSE_CFG.get("predictions_parent_id", ""),
        cluster_keys=cluster_keys_arg(),
    conda:
        "../envs/post_process.yaml"
    resources:
        mem_mb=determine_mem_mb,
    shell:
        """
        python {workflow.basedir}/scripts/manage_synapse_manifest.py \
        --product element_list \
        --parent-id {params.parent_id} \
        --cluster-keys "{params.cluster_keys}" \
        --nested \
        --model-tokens {ELEMENT_LIST_TOKENS} \
        --manifest-out {output} \
        --dry-run {DRY_RUN} \
        --confirm-delete {CONFIRM_DELETE} \
        --confirm-overwrite {CONFIRM_OVERWRITE} \
        {input}
        """


rule manage_candidates_manifest:
    input:
        get_candidates_output_files(),
    output:
        update(os.path.join(MANIFEST_DIR, "candidates_manifest.tsv")),
    params:
        parent_id=SYNAPSE_CFG.get("candidates_parent_id", ""),
        cluster_keys=cluster_keys_arg(),
    conda:
        "../envs/post_process.yaml"
    resources:
        mem_mb=determine_mem_mb,
    shell:
        """
        python {workflow.basedir}/scripts/manage_synapse_manifest.py \
        --product candidates \
        --parent-id {params.parent_id} \
        --cluster-keys "{params.cluster_keys}" \
        --manifest-out {output} \
        --dry-run {DRY_RUN} \
        --confirm-delete {CONFIRM_DELETE} \
        --confirm-overwrite {CONFIRM_OVERWRITE} \
        {input}
        """


rule manage_features_manifest:
    input:
        get_features_output_files(),
    output:
        update(os.path.join(MANIFEST_DIR, "features_manifest.tsv")),
    params:
        parent_id=SYNAPSE_CFG.get("features_parent_id", ""),
        cluster_keys=cluster_keys_arg(),
    conda:
        "../envs/post_process.yaml"
    resources:
        mem_mb=determine_mem_mb,
    shell:
        """
        python {workflow.basedir}/scripts/manage_synapse_manifest.py \
        --product features \
        --parent-id {params.parent_id} \
        --cluster-keys "{params.cluster_keys}" \
        --nested \
        --model-tokens {MODEL_TOKENS} \
        --manifest-out {output} \
        --dry-run {DRY_RUN} \
        --confirm-delete {CONFIRM_DELETE} \
        --confirm-overwrite {CONFIRM_OVERWRITE} \
        {input}
        """
