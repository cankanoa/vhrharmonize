# Installation

## Install with PyPI

1. Create and activate a conda environment.

```bash
conda create -n vhrharmonize -c conda-forge gdal sixs python=3.11
conda activate vhrharmonize
```

2. Install the package with default dependencies.

```bash
pip install vhrharmonize[defaults]
```

3. Install specific dependencies as needed.

```bash
pip install "vhrharmonize[cloud]"
pip install "vhrharmonize[py6s]"
pip install "vhrharmonize[flaash]"
pip install "vhrharmonize[orthorectification]"
pip install "vhrharmonize[pansharpen]"
pip install "vhrharmonize[align]"
pip install "vhrharmonize[spectralmatch]"
pip install "vhrharmonize[docs]"
pip install "vhrharmonize[all]"
```

## Install from source

1. Create and activate a conda environment.

```bash
conda create -n vhrharmonize -c conda-forge gdal sixs python=3.11
conda activate vhrharmonize
```

2. Clone the repository.

```bash
git clone https://github.com/cankanoa/vhrharmonize.git
cd vhrharmonize
```

3. Install the default packages.

```bash
pip install -e '.[defaults]'
```

4. Install specific dependencies as needed.

```bash
pip install -e ".[cloud]"
pip install -e ".[py6s]"
pip install -e ".[flaash]"
pip install -e ".[orthorectification]"
pip install -e ".[pansharpen]"
pip install -e ".[align]"
pip install -e ".[spectralmatch]"
pip install -e ".[docs]"
pip install -e ".[all]"
```

## Py6S and 6S troubleshooting

Install the 6S executable from conda-forge with `conda install conda-forge::sixs`. The `py6s` extra installs the Python interface. If the executable is not available on `PATH`, set `sixs_executable` in the atmospheric correction plugin and run `vhr atmospheric_correction --config recipe.yml`.

5. Verify the entry points if desired.

```bash
vhr workflow --help
vhr fetch_atmosphere --help
vhr atmospheric_correction --help
vhr cloud_mask --help
vhr pansharpen --help
vhr alignment --help
vhr orthorectification --help
vhr global_regression --help
vhr hpc-prepare --help
```

The individual SpectralMatch adapters require SpectralMatch 1.6 or later (below 2.0), including the public `create_footprints`, `postprocess_footprints` and `markov_triangles` functions. Install this version in both the local and HPC environments.
