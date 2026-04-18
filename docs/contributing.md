# Contributing

Contributions are welcome. This document describes how to set up a development environment and the conventions the
project follows.

---

## Development setup

```bash
git clone https://github.com/exterex/clystere.git
cd clystere

# Create a Python virtual environment (Python 3.11+)
python -m venv .venv
source .venv/bin/activate

# Install development tools
pip install pre-commit pytest ruff

# Install and activate the pre-commit hooks
pre-commit install
```

---

## Code style

### Python (`bin/`)

Python scripts are linted and formatted with [Ruff](https://docs.astral.sh/ruff/). The target interpreter is Python
3.11.

### Nextflow (`.nf`, `.config`)

Nextflow files are formatted with [Prettier](https://prettier.io/) via `prettier-plugin-groovy`.

---

## Pre-commit hooks

The following hooks run automatically on `git commit`:

| Hook                  | Purpose                                       |
| --------------------- | --------------------------------------------- |
| `ruff`                | Python linting with auto-fix                  |
| `ruff-format`         | Python formatting                             |
| `prettier`            | Formatting for Nextflow, JSON, YAML, Markdown |
| `trailing-whitespace` | Strip trailing whitespace                     |
| `end-of-file-fixer`   | Ensure files end with a newline               |
| `check-yaml`          | Validate YAML syntax                          |
| `check-json`          | Validate JSON syntax                          |
| `mixed-line-ending`   | Enforce LF line endings                       |

Run all hooks manually:

```bash
pre-commit run --all-files
```

---

## Testing

### Python unit tests

```bash
pytest tests/python/ -v
```

Tests live in `tests/python/` and exercise the pure Python functions in `bin/`. Fixtures are built inline with
`tmp_path` — no external data files needed.

### nf-test (Nextflow module tests)

```bash
# requires Nextflow and Docker
nf-test test --profile docker --verbose
```

Module tests live in `tests/modules/local/<module_name>/main.nf.test` and use a minimal antiSMASH JSON fixture in
`tests/data/`.

### Full pipeline test

```bash
nextflow run . \
    --input assets/samplesheet.csv \
    --outdir results \
    -profile test,docker
```

---

## Branching and pull requests

- `main` is the stable release branch.
- Feature branches: `feature/<description>`
- Bug/hot fixes: `fix/<description>`
- Non-code tasks: `chore/<description>`
- Preparing a release: `release/<description>`
- Open a pull request against `main`; CI must be green before merge.

---

## Releases

Releases follow [Semantic Versioning](https://semver.org/). To cut a release:

1. Update the version in `nextflow.config` (`manifest.version`).
2. Commit with message `release: bump version to vX.Y.Z`.
3. Push a `vX.Y.Z` tag — the `release.yml` workflow creates the GitHub Release automatically.

```bash
git tag v1.1.0
git push origin v1.1.0
```
