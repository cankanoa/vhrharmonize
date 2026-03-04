"""Canonical cloudmask CLI entrypoint."""

import sys

from .cloudmask_on_raster import main

__all__ = ["main"]


if __name__ == "__main__":
    sys.exit(main())
