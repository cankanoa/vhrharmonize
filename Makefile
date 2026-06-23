SHELL := /bin/bash

.PHONY: help clean-package python-build version tag release

help:
	@echo "Python package release helpers:"
	@echo "  make clean-package"
	@echo "  make python-build"
	@echo "  make version version=0.0.2"
	@echo "  make tag version=0.0.2"
	@echo "  make release version=0.0.2"

clean-package:
	rm -rf build dist *.egg-info vhrharmonize.egg-info

python-build: clean-package
	python -m build --sdist --wheel --no-isolation --outdir dist/

tag:
	@if [ -z "$(version)" ]; then \
		echo "Usage: make tag version=0.0.2"; \
		exit 1; \
	fi
	git tag -a v$(version) -m "Version $(version)"
	git push origin v$(version)

version:
	@if [ -z "$(version)" ]; then \
		echo "Usage: make version version=0.0.2"; \
		exit 1; \
	fi
	@echo "Updating pyproject.toml version to $(version)..."
	sed -i.bak "s/^version = .*/version = \"$(version)\"/" pyproject.toml && rm pyproject.toml.bak
	git add pyproject.toml
	git commit -m "Version $(version) released"
	git push origin HEAD

release: version tag
	@echo "Created and pushed release v$(version). GitHub Actions will create the GitHub release and publish to PyPI."
