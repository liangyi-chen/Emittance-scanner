#!/bin/bash
trap "exit 1" INT

# ---------------------------------------------------------------
# Usage:
#   ./emittance_scan.sh --continue   skip points that already have a PNG
#   ./emittance_scan.sh --restart    rerun everything (clears scan.dat)
#   ./emittance_scan.sh              asks which one you want
# ---------------------------------------------------------------
mode=""
case "$1" in
  -c|--continue|--resume) mode="continue" ;;
  -r|--restart|--fresh)   mode="restart" ;;
  "") ;;
  *) echo "Usage: $0 [--continue | --restart]"; exit 1 ;;
esac

if [ -z "$mode" ]; then
  while true; do
    read -r -p "Continue the previous scan (skip finished points) or restart from scratch? [c/r]: " ans
    case "$ans" in
      c|C) mode="continue"; break ;;
      r|R) mode="restart"; break ;;
      *) echo "Please type c or r." ;;
    esac
  done
fi

if [ "$mode" = "restart" ]; then
  echo "Mode: restart (all points will be rerun, scan.dat will be cleared)"
  make clean && make simu save_png || exit 1
  rm -f scan.dat
else
  echo "Mode: continue (points with an existing zy_view PNG will be skipped)"
  make simu save_png || exit 1
fi

# Format a number the same way the C++ programs do for the file names
# (%g, and -0 becomes 0), so "0.3" -> 0.3, "-0.0" -> 0, "100" -> 100
fmt() { awk -v v="$1" 'BEGIN { x = v + 0; if (x == 0) x = 0; printf "%g", x }'; }

for y in $(seq -1.6 0.2 1.6); do
  for V in $(seq -1000 200 1000); do
    ty=$(fmt "$y")
    tv=$(fmt "$V")
    png="zy_view_${ty}_${tv}.png"

    if [ "$mode" = "continue" ]; then
      if [ -s "$png" ]; then
        echo "y = $ty mm, V = $tv V: Data found. Skip this step."
        continue
      fi
    fi

    ./simu.run 4 $y $V
    ./save_png.run $y $V
  done
done
