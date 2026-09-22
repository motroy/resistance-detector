#!/usr/bin/env python3
"""Compatibility shim: the builder now lives in the installable package.

Run ``python3 -m fos_cazavi.build_data`` directly, or keep using this path.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from fos_cazavi.build_data import main  # noqa: E402

if __name__ == '__main__':
    main()
