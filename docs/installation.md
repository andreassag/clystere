# Installation

## Requirements

| Dependency                           | Minimum version | Notes                                                   |
| ------------------------------------ | --------------- | ------------------------------------------------------- |
| [Nextflow](https://www.nextflow.io/) | 23.04.0         | `NXF_VER` env var can pin the version                   |
| Java                                 | 11              | Required by Nextflow                                    |
| Container engine                     | —               | One of Docker, Singularity, Apptainer, Podman, or Conda |

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

## 3. antiSMASH database

The pipeline can download the database automatically on first run. To pre-download it manually:

<!-- markdownlint-disable MD046 -->

=== "Docker"

    ```bash
    # pull the antiSMASH container and run the download helper
    docker pull antismash/standalone:8.0.4
    docker run --rm -v /path/to/db:/db \
        antismash/standalone:8.0.4 \
        antismash-download-databases /db
    ```

=== "Conda"

    ```bash
    # create an antiSMASH environment and run the download helper
    conda create -n antismash -c conda-forge -c bioconda antismash=8
    conda activate antismash
    antismash-download-databases
    ```

<!-- markdownlint-enable MD046 -->

Pass the database path to the pipeline with `--antismash_db /path/to/db`. If the path is absent or empty, the pipeline
downloads the database there automatically before running antiSMASH.

deepBGC model data is also downloaded automatically on first run unless `--deepbgc_data_dir` points to an existing
directory. Ensure outbound network access for the initial deepBGC run.

---

## 4. Clone the repository (optional)

Running `nextflow run andreassag/clystere` pulls the pipeline automatically from GitHub. To work with a local copy:

```bash
git clone https://github.com/andreassag/clystere.git
cd clystere
nextflow run . --input samplesheet.csv --outdir results -profile docker
```
