#!/bin/bash
# Assemble frames from scripts/web_frames.jl into a looping animated WebP (needs img2webp from libwebp).
# Each frame shows for 100 ms; the last frame holds for 3 s before the loop restarts.
#
# Usage: scripts/make_webp.sh <mode> [quality]     e.g. scripts/make_webp.sh white 75
#   reads output/web_frames/<mode>/frame_*.png, writes output/logo_<mode>.webp
set -euo pipefail
mode=$1; q=${2:-75}
cd "$(dirname "$0")/.."
files=(output/web_frames/$mode/frame_*.png)
n=${#files[@]}
args=()
for ((i = 0; i < n; i++)); do
  if (( i == n - 1 )); then args+=(-d 3000 "${files[$i]}"); else args+=(-d 100 "${files[$i]}"); fi
done
img2webp -loop 0 -lossy -q "$q" -m 4 "${args[@]}" -o "output/logo_$mode.webp"
echo "Saved output/logo_$mode.webp"
