"""Canonical MODIS water vapor CLI entrypoint."""

import sys

from .fetch_modis_water_vapor_gee import main

__all__ = ["main"]


if __name__ == "__main__":
    sys.exit(main())
