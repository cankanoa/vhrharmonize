SHELL := /bin/bash

.PHONY: help clean-package python-build build-qgis upload-qgis python-version qgis-version tag release

help:
	@echo "Python package release helpers:"
	@echo "  make clean-package"
	@echo "  make python-build"
	@echo "  make build-qgis"
	@echo "  make upload-qgis"
	@echo "  make python-version version=0.0.2"
	@echo "  make qgis-version version=0.0.2"
	@echo "  make tag version=0.0.2"
	@echo "  make release version=0.0.2"

clean-package:
	rm -rf build dist *.egg-info vhrharmonize.egg-info

python-build: clean-package
	python -m build --sdist --wheel --no-isolation --outdir dist/

build-qgis:
	python qgis_vhrharomize/build_plugin.py

upload-qgis: build-qgis
	python qgis_vhrharomize/plugin_upload.py qgis_vhrharomize.zip

tag:
	@if [ -z "$(version)" ]; then \
		echo "Usage: make tag version=0.0.2"; \
		exit 1; \
	fi
	git tag -a v$(version) -m "Version $(version)"
	git push origin v$(version)

python-version:
	@if [ -z "$(version)" ]; then \
		echo "Usage: make python-version version=0.0.2"; \
		exit 1; \
	fi
	@echo "Updating pyproject.toml version to $(version)..."
	sed -i.bak "s/^version = .*/version = \"$(version)\"/" pyproject.toml && rm pyproject.toml.bak
	git add pyproject.toml
	git commit -m "Update Python version to $(version)"
	git push origin HEAD

qgis-version:
	@if [ -z "$(version)" ]; then \
		echo "Usage: make qgis-version version=0.0.2"; \
		exit 1; \
	fi
	@echo "Updating QGIS metadata version to $(version)..."
	sed -i.bak "s/^version=.*/version=$(version)/" qgis_vhrharomize/metadata.txt && rm qgis_vhrharomize/metadata.txt.bak
	git add qgis_vhrharomize/metadata.txt
	git commit -m "Update QGIS version to $(version)"
	git push origin HEAD

release: python-version qgis-version tag
	@echo "Created and pushed release v$(version). GitHub Actions will create the GitHub release and publish to PyPI."
