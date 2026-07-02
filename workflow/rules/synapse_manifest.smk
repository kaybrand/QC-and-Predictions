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
"""

SYNAPSE_CFG = config.get("synapse", {})
DRY_RUN = str(SYNAPSE_CFG.get("dry_run", True))
CONFIRM_DELETE = str(SYNAPSE_CFG.get("confirm_delete", False))
CONFIRM_OVERWRITE = str(SYNAPSE_CFG.get("confirm_overwrite", False))
MANIFEST_DIR = os.path.join(RESULTS_DIR_BASE, "synapse_manifests")
MODEL_TOKENS = "multiome_powerlaw_v3,scATAC_powerlaw_v3"


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


rule manage_filtered_data_manifest:
    input:
        get_filtered_data_files(),
    output:
        os.path.join(MANIFEST_DIR, "filtered_data_manifest.tsv"),
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
        get_reformat_output_files(),
    output:
        os.path.join(MANIFEST_DIR, "predictions_manifest.tsv"),
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


rule manage_candidates_manifest:
    input:
        get_candidates_output_files(),
    output:
        os.path.join(MANIFEST_DIR, "candidates_manifest.tsv"),
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
        os.path.join(MANIFEST_DIR, "features_manifest.tsv"),
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
