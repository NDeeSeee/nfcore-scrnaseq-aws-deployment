#!/bin/bash
# =============================================================================
# install_tools_bone10x.sh — run ONCE on a login node (needs outbound network).
# Installs simpleaf + piscem + alevin-fry into a micromamba env under
# $WD/tools_bone10x and seeds a fresh ALEVIN_FRY_HOME. Replaces the March tools
# that lived in the deactivated pavb5f account. Touches nothing that exists.
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/bone10x_refseq_env.sh"

mkdir -p "$TOOLROOT/bin"
MM="$TOOLROOT/bin/micromamba"
if [[ ! -x "$MM" ]]; then
    echo "== fetching micromamba =="
    curl -fsSL https://micro.mamba.pm/api/micromamba/linux-64/latest \
        | tar -xj -C "$TOOLROOT" bin/micromamba
fi

if [[ ! -x "$TOOLBIN/simpleaf" ]]; then
    echo "== creating env (alevin-fry pinned to the March version) =="
    "$MM" create -y -r "$TOOLROOT" -n af \
        --strict-channel-priority -c conda-forge -c bioconda \
        simpleaf piscem alevin-fry=0.11.2 'salmon>=1.10,<2'
    # salmon is never used (piscem maps), but simpleaf 0.20 version-checks any
    # salmon on PATH and rejects 2.x, which bioconda otherwise pulls in.
fi

echo "== simpleaf home: $ALEVIN_FRY_HOME =="
mkdir -p "$ALEVIN_FRY_HOME"
# Reuse the permit lists March already downloaded so compute nodes stay offline.
OLD_HOME="$WD/.alevin_fry_home"
[[ -d "$ALEVIN_FRY_HOME/plist" ]]            || cp -r "$OLD_HOME/plist" "$ALEVIN_FRY_HOME/"
[[ -f "$ALEVIN_FRY_HOME/chemistries.json" ]] || cp "$OLD_HOME/chemistries.json" "$ALEVIN_FRY_HOME/"
simpleaf set-paths

for t in simpleaf piscem alevin-fry; do
    printf '%-11s %s\n' "$t" "$("$t" --version 2>&1 | head -1)"
done
cat "$ALEVIN_FRY_HOME/simpleaf_info.json"
echo "Install OK."
