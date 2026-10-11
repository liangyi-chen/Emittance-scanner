# Emittance scan with Allison-type scanner in IBSimu

This IBSimu project simulates emittance measurement of an ion beam with Allison-type emittance scanner. The goal is to establish a first-principle simulation for the emittance scanner. With the inclusion of magnetic field, the final goal is to quantify the effect of external magnetic field on measured emittance.

## Use of simu.cpp

simu.cpp is used to run simulation. It takes three arguments: number of CPU cores (default is 4), y position of emittance scanner (unit is mm), and deflection voltage of plate electrode (unit is V). After compiling, run the simulation with:

`./simu.run 4 0.2 200`

Data files are saved in `geom/epot/pdb_yposition_voltage.dat` format. The current on the Faraday cup of the emittance scanner is also appended to `scan.dat`.

## Use of plot.cpp

plot.cpp plots the result. It takes the filenames of the simulation data files as arguments. A typical command to run plot.cpp is:

`./plot.run geom_0.2_200.dat epot_0.2_200.dat pdb_0.2_200.dat`

## Use of save_png.cpp

save_png.cpp skips the gtk plotter window and directly prints a png image showing the yz view of the beam through the scanner. The syntax to run save_png.cpp is 

`./save_pgn.run 0.2 200`

## Automatic scan sequence using emittance_scan.sh

Use the shell script to run through scanning the y position and deflection voltage of the emittance scanner. The scans are in range `(seq start step stop)` for both parameters defined in the shell script.
