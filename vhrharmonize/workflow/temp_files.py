"""Redirect explicitly selected file arguments for one plugin invocation."""

from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
import shutil
from uuid import uuid4


class TempFileSwapError(OSError):
    """A core file-transfer failure; never an attributable image failure."""


@dataclass
class TempFileSwap:
    names: tuple[str, ...]
    directory: str
    input_names: frozenset[str]
    output_names: frozenset[str]
    requirements: tuple[str, ...] = ()

    @contextmanager
    def apply(self, params):
        """Keep original destinations outside the call; copy outputs back on success."""
        files, retained = {}, set()

        def redirect(value, *, is_input, is_output):
            if value is None:
                return None
            if isinstance(value, list):
                return [redirect(item, is_input=is_input, is_output=is_output) for item in value]
            if not isinstance(value, str) or not value:
                raise TempFileSwapError("core:swap_for_temp_file accepts file paths or lists of file paths")
            if value not in files:
                source = Path(value)
                folder = Path(self.directory) / uuid4().hex
                local = folder / source.name
                try:
                    folder.mkdir(parents=True)
                    files[value] = (local, is_output)
                    if is_input or source.exists():
                        if not source.is_file():
                            raise TempFileSwapError(f"core:swap_for_temp_file requires a regular file: {value}")
                        shutil.copy2(source, local)
                    if is_input:
                        for filename in self.requirements:
                            companion = Path(filename)
                            if (companion != source and companion.parent == source.parent
                                    and companion.name.startswith(source.stem + ".")):
                                shutil.copy2(companion, folder / companion.name)
                except OSError as exc:
                    raise TempFileSwapError(f"Could not prepare temporary file {local}: {exc}") from exc
            else:
                local, output = files[value]
                files[value] = (local, output or is_output)
            return str(files[value][0])

        try:
            swapped = dict(params)
            for name in self.names:
                if name not in params:
                    raise TempFileSwapError(f"core:swap_for_temp_file parameter is not supplied: {name}")
                swapped[name] = redirect(params[name], is_input=name in self.input_names,
                                         is_output=name in self.output_names)
            originals = {str(local): original for original, (local, _) in files.items()}
            yield swapped, originals
            outputs = [(original, local) for original, (local, output) in files.items() if output]
            for original, local in outputs:
                if not local.is_file():
                    raise TempFileSwapError(f"Temporary output was not created for {original}: {local}")
            retained = {local.parent for _, local in outputs}
            for original, local in outputs:
                try:
                    Path(original).parent.mkdir(parents=True, exist_ok=True)
                    # NAS destinations may allow writes but reject timestamp changes.
                    shutil.copyfile(local, original)
                except OSError as exc:
                    raise TempFileSwapError(f"Could not copy temporary output to {original}; retained at {local}") from exc
                retained.remove(local.parent)
        finally:
            for local, _ in files.values():
                if local.parent not in retained:
                    try:
                        shutil.rmtree(local.parent)
                    except OSError as exc:
                        raise TempFileSwapError(f"Could not remove temporary folder {local.parent}: {exc}") from exc
