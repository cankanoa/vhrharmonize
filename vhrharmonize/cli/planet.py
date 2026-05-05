"""Planet CLI scaffold."""

import argparse
import sys
from typing import Optional


def main(argv: Optional[list[str]] = None) -> int:
    """Run the placeholder Planet CLI.
    Args:
        argv: Optional command line arguments.
    Returns:
        Process exit code.
    """
    parser = argparse.ArgumentParser(
        description="Planet workflow support is not implemented yet."
    )
    parser.parse_args(argv)
    parser.exit(
        status=1,
        message="Planet workflow support is not implemented yet.\n",
    )


if __name__ == "__main__":
    sys.exit(main())
