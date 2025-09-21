#!/usr/bin/env bash
set -euo pipefail

ENV_FILE="$(dirname "$0")/.env"

set -a
. "$ENV_FILE"
set +a

WSL_PROJECT_ROOT="$(wslpath -a "$WIN_PROJECT_ROOT" 2>/dev/null || echo "$WIN_PROJECT_ROOT")"
export WSL_PROJECT_ROOT
export WSL_REACTIONS="$WSL_PROJECT_ROOT/$REACTIONS_DIR"
export WSL_REACTION_DIR="$WSL_REACTIONS/$REACTION"
export WSL_INPUTS="$WSL_REACTION_DIR/$INPUTS_SUBDIR"
export WSL_RUNS="$WSL_REACTION_DIR/$RUNS_SUBDIR"