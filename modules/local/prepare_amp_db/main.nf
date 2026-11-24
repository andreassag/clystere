/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PREPARE_AMP_DB
    Downloads or indexes the AMPcombi2 reference database once per run so parallel
    AMPCOMBI2_PARSETABLES tasks never race on mmseqs createdb or duplicate downloads.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process PREPARE_AMP_DB {
    tag "${opt_amp_db}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/06/063288a19d4f267756c136b1ac8bd06ef89ecee62bf6fcc5b27732848f118b01/data'
        : 'community.wave.seqera.io/library/ampcombi:3.0.0--3b8b8b73f38675a2'}"

    input:
    val opt_amp_db
    path existing_db_dir, stageAs: 'existing_db'
    val trigger

    output:
    path "amp_${opt_amp_db}_database", emit: db
    tuple val("${task.process}"), val('ampcombi'), eval("ampcombi --version | sed 's/ampcombi //'"), emit: versions_prepare_amp_db, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def has_existing = existing_db_dir ? 'true' : 'false'
    """
    mkdir -p "amp_${opt_amp_db}_database"
    if [[ "${has_existing}" == "true" ]]; then
        cp -rL "${existing_db_dir}/." "amp_${opt_amp_db}_database/"
    fi

    python3 - "amp_${opt_amp_db}_database" "${opt_amp_db}" "${has_existing}" "${task.cpus}" <<'PY'
    import shutil
    import sys
    from pathlib import Path
    import ampcombi

    pkg_dir = str(Path(ampcombi.__file__).parent)
    if pkg_dir not in sys.path:
        sys.path.insert(0, pkg_dir)

    import amp_database

    db_dir = sys.argv[1]
    database = sys.argv[2]
    has_existing = sys.argv[3] == "true"
    threads = int(sys.argv[4])

    if not has_existing:
        amp_database.download_ref_db(db_dir, database, threads)

    db_path = Path(db_dir)
    mmseqs_dir = db_path / "mmseqs2"
    ref_db_dbtype = mmseqs_dir / "ref_DB.dbtype"
    if mmseqs_dir.exists() and not ref_db_dbtype.is_file():
        shutil.rmtree(mmseqs_dir)

    if not ref_db_dbtype.is_file():
        fasta_files = list(db_path.glob("*.fasta"))
        if not fasta_files:
            fa_files = sorted(db_path.glob("*.fa"))
            if not fa_files:
                raise RuntimeError(
                    f"AMP reference database directory '{db_dir}' contains no .fasta/.fa files and no mmseqs2/ref_DB."
                )
            for fa in fa_files:
                link = fa.with_suffix(".fasta")
                if not link.exists():
                    link.symlink_to(fa.name)
        amp_database.create_mmseqs_ref_db(db_dir)

    if not ref_db_dbtype.is_file():
        raise RuntimeError(f"Failed to prepare MMseqs2 reference database in '{db_dir}'.")
    PY
    """

    stub:
    """
    mkdir -p "amp_${opt_amp_db}_database/mmseqs2"
    touch "amp_${opt_amp_db}_database/general_amps.fasta"
    touch "amp_${opt_amp_db}_database/general_amps.txt"
    touch "amp_${opt_amp_db}_database/mmseqs2/ref_DB"
    touch "amp_${opt_amp_db}_database/mmseqs2/ref_DB.dbtype"
    touch "amp_${opt_amp_db}_database/mmseqs2/ref_DB_h"
    touch "amp_${opt_amp_db}_database/mmseqs2/ref_DB_h.dbtype"
    touch "amp_${opt_amp_db}_database/mmseqs2/ref_DB_h.index"
    touch "amp_${opt_amp_db}_database/mmseqs2/ref_DB.index"
    touch "amp_${opt_amp_db}_database/mmseqs2/ref_DB.lookup"
    touch "amp_${opt_amp_db}_database/mmseqs2/ref_DB.source"
    """
}
