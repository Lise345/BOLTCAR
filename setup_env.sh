#!/usr/bin/env bash
set -euo pipefail

VENV_DIR="${VENV_DIR:-.venv}"
PYTHON_BIN="${PYTHON:-python3}"

# (Optional) if you know you’ll need a compiler for wheels with C extensions, uncomment:
# module load foss/2023a

echo "→ Creating virtual environment in ${VENV_DIR}"
"${PYTHON_BIN}" -m venv "${VENV_DIR}"

# shellcheck disable=SC1090
source "${VENV_DIR}/bin/activate"

echo "→ Upgrading pip tooling"
python -m pip install --upgrade pip setuptools wheel

echo "→ Installing requirements"
pip install -r requirements.txt

echo "→ Freezing a lock file (optional but nice for reproducibility)"
pip freeze --local > requirements.lock.txt

echo "✅ Environment ready at ${VENV_DIR}"
