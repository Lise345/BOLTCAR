#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
VENV_DIR="${VENV_DIR:-$SCRIPT_DIR/.venv}"
if [ ! -x "$VENV_DIR/bin/python" ]; then
  python3 -m venv "$VENV_DIR"
  "$VENV_DIR/bin/pip" install --upgrade pip setuptools wheel
  [ -f "$SCRIPT_DIR/requirements.txt" ] && "$VENV_DIR/bin/pip" install -r "$SCRIPT_DIR/requirements.txt"
fi
# shellcheck disable=SC1090
source "$VENV_DIR/bin/activate"
