"""Read snapshot comment headers followed by an ordinary CSV table."""

import csv
from itertools import chain
import math


def read_snapshot_csv(path):
    """Return (metadata, column_names, rows); table values remain strings.

    New snapshots start with '# step=...', '# time=...', '# dt=...'. Old
    files without a comment header return an empty metadata dictionary.
    """
    metadata = {}
    with path.open(newline="") as stream:
        for line in stream:
            if not line.strip():
                continue
            if not line.startswith("#"):
                reader = csv.DictReader(chain((line,), stream))
                fields = reader.fieldnames or []
                rows = list(reader)
                break
            key, separator, value = line[1:].strip().partition("=")
            key = key.strip()
            if key in ("step", "time", "dt"):
                if not separator or key in metadata:
                    raise ValueError(f"Invalid or duplicate metadata '{key}' in {path}")
                metadata[key] = int(value) if key == "step" else float(value)
        else:
            raise ValueError(f"Missing CSV table in {path}")
    if metadata:
        if set(metadata) != {"step", "time", "dt"}:
            raise ValueError(f"Incomplete step/time/dt header in {path}")
        if metadata["step"] < 0 or any(not math.isfinite(metadata[key]) or metadata[key] < 0
                                       for key in ("time", "dt")):
            raise ValueError(f"Invalid snapshot step/time/dt in {path}")
    return metadata, fields, rows
