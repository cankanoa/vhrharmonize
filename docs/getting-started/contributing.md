# Contributing

## Structure

Basic project layout:

- `vhrharmonize/cli`: packaged command-line entry points
- `vhrharmonize/plugins`: self-contained processing implementations and plugin adapters
- `vhrharmonize/workflow`: planning, paths, metadata and plugin registration
- `vhrharmonize/io`: shared raster and workflow IO helpers
- `configs`: example configuration files
- `docs`: MkDocs source
- `tests`: tests

## Build Docs Locally

```bash
pip install -e ".[docs]"
mkdocs serve -f docs/mkdocs.yml
```

```bash
mkdocs build -f docs/mkdocs.yml
```

## Releases

```bash
make release version=0.0.2
```

This updates `pyproject.toml`, commits the version bump, creates the `v0.0.2` tag, and pushes it. GitHub Actions then creates the GitHub release and publishes to PyPI.

## How to Contribute

1. Open an issue first if the change affects behavior, API, or workflow design.
2. Fork the repository.
3. Create a branch in your fork.
4. Make the change.
5. Open a pull request from your fork to the main repository.

CLI wrappers are generated from public function signatures/docstrings by `vhrharmonize.cli.functions`. Keep validation and defaults in the Python APIs, and add a CLI/API parity test when changing behavior.

See [Adding plugins](adding-plugins.md) for registration and the plugin CLI.
