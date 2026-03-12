#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
yaml_path="${1:-$repo_root/.secrets/openai.yaml}"

if [[ ! -f "$yaml_path" ]]; then
  echo "Missing YAML file: $yaml_path" >&2
  echo "Create it from .secrets/openai.yaml.example first." >&2
  exit 1
fi

key="$(
python3 - "$yaml_path" <<'PY'
import sys

try:
    import yaml
except ImportError:
    raise SystemExit("PyYAML is required. Install with: python3 -m pip install --user pyyaml")

path = sys.argv[1]
with open(path, "r", encoding="utf-8") as f:
    data = yaml.safe_load(f) or {}

try:
    key = data["openai"]["api_key"]
except Exception:
    raise SystemExit(f"Could not find openai.api_key in {path}")

if not isinstance(key, str) or not key.strip():
    raise SystemExit(f"openai.api_key is empty in {path}")

print(key.strip())
PY
)"

export OPENAI_API_KEY="$key"
echo "OPENAI_API_KEY loaded from $yaml_path"
