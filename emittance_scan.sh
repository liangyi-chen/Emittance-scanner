#!/bin/bash
trap "exit 1" INT
make clean && make simu save_png || exit 1
rm -f scan.dat
for y in $(seq -0.5 0.1 0.5); do
  for V in $(seq -300 100 300); do
    ./simu.run 4 $y $V
    ./save_png.run $y $V
  done
done
