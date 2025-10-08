#!/usr/bin/env bash
set -euo pipefail

# Resolve this script's directory (works from anywhere)
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

# Where to keep the venv (override with: export VENV_DIR=/path/to/.venv)
VENV_DIR="${VENV_DIR:-$SCRIPT_DIR/.venv}"
VENV_BIN="$VENV_DIR/bin"

# Bootstrap venv if missing
if [ ! -x "$VENV_BIN/python" ]; then
  echo "→ Creating venv at $VENV_DIR"
  python3 -m venv "$VENV_DIR"
  "$VENV_BIN/pip" install --upgrade pip setuptools wheel
  if [ -f "$SCRIPT_DIR/requirements.txt" ]; then
    echo "→ Installing requirements.txt"
    "$VENV_BIN/pip" install -r "$SCRIPT_DIR/requirements.txt"
  fi
fi

# Activate venv for tools that expect activation
# shellcheck disable=SC1090
source "$VENV_BIN/activate"

# Helpful warning if Gaussian isn't on PATH (optional)
if ! command -v g16 >/dev/null 2>&1 && ! command -v g09 >/dev/null 2>&1; then
  echo "⚠️  Gaussian not found in PATH. Load its module if you need it." >&2
fi

# If no args: just drop into a shell with venv active
if [ "$#" -eq 0 ]; then
  echo "Activated $VENV_DIR. No command given; opening a shell."
  exec "${SHELL:-/bin/bash}"
fi

# If first arg isn't a command but looks like a Python script, or a file,
# auto-prepend 'python'
first="$1"
if ! command -v "$first" >/dev/null 2>&1; then
  if [[ "$first" == *.py ]] || [ -f "$first" ]; then
    set -- python "$@"
  fi
fi

exec "$@"
