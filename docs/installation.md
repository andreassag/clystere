# Installation

## Requirements

| Dependency                           | Minimum version | Notes                                                   |
| ------------------------------------ | --------------- | ------------------------------------------------------- |
| [Nextflow](https://www.nextflow.io/) | 25.10.0         | `NXF_VER` env var can pin the version                   |
| Java                                 | 17              | Required by Nextflow (Java 17 or higher)                |
| Container engine                     | —               | One of Docker, Apptainer, Singularity, Podman, or Conda |

---

## 1. Install Nextflow

```bash
# via the official installer
curl -s https://get.nextflow.io | bash
mv nextflow ~/.local/bin/

# verify
nextflow -version
```

Or with Conda/Mamba:

```bash
mamba create -n nf nextflow
conda activate nf
```

---

## 2. Install a container engine

<!-- markdownlint-disable MD046 -->

=== "Docker"

    Follow the [Docker installation guide](https://docs.docker.com/get-docker/).
    Verify with `docker --version`.

=== "Apptainer"

    Follow the [Apptainer installation guide](https://apptainer.org/docs/user/latest/quick_start.html).
    Verify with `apptainer --version`.

=== "Conda"

    Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or Mamba.
    All tool environments are managed per-process by Nextflow.

<!-- markdownlint-enable MD046 -->

---

## 3. Reference databases

Reference databases are downloaded and cached automatically under `--clystere_db_dir`
(`${launchDir}/.clystere/databases`) on first run when explicit database paths are not provided:

- **antiSMASH** (`--antismash_db`): downloaded via `ANTISMASH_ANTISMASHDOWNLOADDATABASES` into `antismash_db/`.
- **deepBGC** (`--deepbgc_data_dir`): downloaded via `DEEPBGC_DOWNLOAD` into `deepbgc_db/` when `--deepbgc_run` is
  enabled.
- **BiG-SLiCE** (`--bigslice_db`): downloaded via `BIGSLICE_DOWNLOADDB` into `bigslice-models/` when `--bigslice_run` is
  enabled.

To pre-download the antiSMASH database manually:

<!-- markdownlint-disable MD046 -->

=== "Docker"

    ```bash
    docker pull quay.io/nf-core/antismash:8.0.1--pyhdfd78af_0
    docker run --rm -v /path/to/db:/db \
        quay.io/nf-core/antismash:8.0.1--pyhdfd78af_0 \
        download-antismash-databases --database-dir /db
    ```

=== "Conda"

    ```bash
    conda create -n antismash -c conda-forge -c bioconda antismash=8.0.1
    conda activate antismash
    download-antismash-databases --database-dir /path/to/db
    ```

<!-- markdownlint-enable MD046 -->

---

## 4. Clone the repository (optional)

Running `nextflow run andreassag/clystere` pulls the pipeline automatically from GitHub. To work with a local copy:

```bash
git clone https://github.com/andreassag/clystere.git
cd clystere
nextflow run . --input assets/samplesheet.csv --outdir results -profile docker
```
