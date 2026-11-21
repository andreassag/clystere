# Contributing

Contributions are welcome. This document describes how to set up a development environment and the conventions the
project follows.

---

## Development setup

```bash
git clone https://github.com/andreassag/clystere.git
cd clystere

# Create a Python virtual environment (Python 3.11+)
python -m venv .venv
source .venv/bin/activate

# Install development tools
pip install pre-commit pytest ruff biopython pandas jsonschema
npm install

# Install pre-commit git hooks
pre-commit install
```

---

## Code style

### Python (`bin/`)

Python scripts are linted and formatted with [Ruff](https://docs.astral.sh/ruff/). The target interpreter is Python
3.11.

```bash
ruff check .
ruff format --check .
```

### Nextflow (`.nf`, `.config`) & Documentation

Nextflow files are validated with `nextflow lint .`, while Markdown, YAML, JSON, and file conventions follow
EditorConfig and Prettier standards:

```bash
nextflow lint .
npx editorconfig-checker
npx prettier --check .
```

---

## Pre-commit hooks

Standard [pre-commit](https://pre-commit.com/) hooks defined in `.pre-commit-config.yaml` run automatically on
`git commit`:

| Check               | Tool               | Purpose                          |
| :------------------ | :----------------- | :------------------------------- |
| Python Linting      | `ruff`             | Code style and import order      |
| Python Formatting   | `ruff-format`      | Code layout consistency          |
| Document Formatting | `prettier`         | Markdown, JSON, and YAML styling |
| File Hygiene        | `pre-commit-hooks` | Whitespace, EOF, YAML/JSON check |

Run all hooks and repository checks manually:

```bash
pre-commit run --all-files
npm run check:all
```

---

## Testing

### Python unit tests

```bash
pytest tests/python/ -v
```

Tests live in `tests/python/` and exercise the Python utilities in `bin/`.

### nf-test (Nextflow pipeline, subworkflow & module tests)

```bash
# Run all nf-test suites (integration, configs, subworkflows, and local modules)
nf-test test tests/ --profile apptainer --verbose

# Run specific suites
nf-test test tests/modules/local/summary --profile apptainer --verbose
nf-test test tests/modules/local/combgc --profile apptainer --verbose
nf-test test tests/modules/local/bigscape --profile apptainer --verbose
nf-test test tests/subworkflows/ --profile apptainer --verbose
```

### Full pipeline test

```bash
# Fast compact BGC test profile
nextflow run . -profile test,apptainer

# Full-genome benchmark with BiG-SCAPE
nextflow run . -profile test_full,apptainer
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

1. Update the version in `nextflow.config` (`manifest.version`) and `package.json`.
2. Commit with message `release: bump version to X.Y.Z`.
3. Push a `X.Y.Z` tag — the `release.yml` workflow creates the GitHub Release automatically.

```bash
git tag 1.1.0
git push origin 1.1.0
```
