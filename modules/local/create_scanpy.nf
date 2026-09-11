process CREATE_SCANPY {
    tag "${meta.id}"
    publishDir "${params.outdir}/${meta.id}/scanpy", mode: params.publish_mode

    container params.scanpy_container

    input:
    tuple val(meta), path(inputs)

    output:
    tuple val(meta), path("${meta.id}.scanpy.h5ad"), emit: object

    script:
    """
    # Apptainer/Singularity may expose HOME and the image filesystem read-only.
    # Keep Python/Matplotlib/Numba caches inside this disposable task directory.
    export XDG_CACHE_HOME="\$PWD/.cache"
    export MPLCONFIGDIR="\$PWD/.cache/matplotlib"
    export NUMBA_CACHE_DIR="\$PWD/.cache/numba"
    mkdir -p "\$MPLCONFIGDIR" "\$NUMBA_CACHE_DIR"

    vdj_args=()
    if [[ -f vdj_calls.tsv && -f vdj_receptors.tsv ]]; then
        vdj_args+=(--vdj-calls vdj_calls.tsv --vdj-receptors vdj_receptors.tsv)
    fi

    python3 \$(command -v create_scanpy.py) \
        --exonic exonic \
        --intronic intronic \
        --project '${meta.id}' \
        --output '${meta.id}.scanpy.h5ad' \
        "\${vdj_args[@]}"
    """
    stub:
    """
    touch '${meta.id}.scanpy.h5ad'
    """

}
