"""Canonical pansharpen CLI entrypoint."""

import sys

from .pansharpen_from_orthos import main

__all__ = ["main"]


if __name__ == "__main__":
    sys.exit(main())
