# Test Layout

- `tests/unit/`: fast unit tests by package layer
- `tests/integration/cli/`: CLI behavior and argument/validation tests
- `tests/integration/workflows/`: end-to-end workflow tests (dataset fixtures)

Suggested next additions:

1. provider discovery/metadata parsing tests
2. pipeline config load/merge tests
3. atmospheric adapter contract tests (`FLAASHCorrector`, `Py6SCorrector`)
4. normalization adapter contract tests (`spectralmatch_rrn`)
