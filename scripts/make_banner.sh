#!/bin/bash
# Build the web banner files from a finished simulation run:
#   output/banner/trails.webm       trails-only video, VP9 with alpha (Chrome, Firefox, Edge)
#   output/banner/trails_hevc.mov   trails-only video, HEVC with alpha (Safari; needs macOS VideoToolbox)
#   output/banner/letters.webp      lettering, 2200 px wide, transparent
#   output/banner/still.webp        one trails frame at full development (reduced-motion fallback)
#   output/banner/box_curve.json    background-box darkness per video frame
# Optionally copies them into DEST (e.g. the website's img/banner/).
#
# Usage: scripts/make_banner.sh [DEST]
#   env: RUN (default heated_text_694x231x8), PYTHON (default python3), N (particles, default 5000)
set -euo pipefail
cd "$(dirname "$0")/.."
RUN=${RUN:-heated_text_694x231x8}; PY=${PYTHON:-python3}; N=${N:-5000}
OUT=output/banner; FR=output/web_frames/banner_trails
mkdir -p "$OUT"

echo "1/5 exporting surface velocities (days 0–28)"
julia --project=. scripts/export_velocity.jl "$RUN" surface 0 28 1100

echo "2/5 rendering particle trails"
rm -rf "$FR"
"$PY" scripts/particle_trails.py output/velocity_data_surface "$FR" \
    --n "$N" --trail 48 --life 96,384 --taper 1 --fade 24 --soft --vmin 0.05 --lw 0.9 --alpha 0.6 \
    --color white --edgefade 0.08,0.2 --rampup 10 --stopseed 20 --noletters

echo "3/5 encoding videos"
ffmpeg -v error -y -framerate 10 -i "$FR/frame_%04d.png" -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2:color=0x00000000" \
    -c:v libvpx-vp9 -pix_fmt yuva420p -b:v 0 -crf 45 -row-mt 1 -deadline good -cpu-used 1 -g 240 "$OUT/trails.webm"
ffmpeg -v error -y -framerate 10 -i "$FR/frame_%04d.png" -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2:color=0x00000000,format=bgra" \
    -c:v hevc_videotoolbox -alpha_quality 0.5 -q:v 40 -tag:v hvc1 -movflags +faststart "$OUT/trails_hevc.mov"

echo "4/5 lettering, still frame, box curve"
julia --project=. scripts/export_letters.jl "$RUN" 2200 "$OUT/letters.png"
cwebp -quiet -lossless "$OUT/letters.png" -o "$OUT/letters.webp" && rm "$OUT/letters.png"
cwebp -quiet -lossless "$FR/frame_0152.png" -o "$OUT/still.webp"            # day 19, just before thinning starts
"$PY" scripts/box_curve.py "$FR" "$OUT/box_curve.json"

echo "5/5 done:"; ls -la "$OUT"
if [ $# -ge 1 ]; then mkdir -p "$1"; cp "$OUT"/{trails.webm,trails_hevc.mov,letters.webp,still.webp,box_curve.json} "$1"/ && echo "copied to $1"; fi
