# Emittance scan with Allison-type scanner in IBSimu

This IBSimu project simulates emittance measurement of an ion beam with Allison-type emittance scanner. The goal is to establish a first-principle simulation for the emittance scanner. With the inclusion of magnetic field, the final goal is to quantify the effect of external magnetic field on measured emittance.

## Use of simu.cpp

simu.cpp is used to run simulation. It takes three arguments: number of CPU cores (default is 4), y position of emittance scanner (unit is mm), and deflection voltage of plate electrode (unit is V). After compiling, run the simulation with:

`./simu.run 4 0.2 200`

Data files are saved in `geom/epot/pdb_yposition_voltage.dat` format. The current on the Faraday cup of the emittance scanner is also appended to `scan.dat`.

### Simulation parameters that are defined within simu.cpp

Set the z position of the emittance scanner with

```c++
const double scanner_z_position = x.0e-3;
```
Set the magnetic field scaling factor inside the emittance scanner with

```c++
double bfield_inside_factor = 1.0;
```

Set the space charge compensation in the free drift space with

```c++
double sc_comp = 0.2;
```

## Use of plot.cpp

plot.cpp plots the result. It takes the filenames of the simulation data files as arguments. A typical command to run plot.cpp is:

`./plot.run geom_0.2_200.dat epot_0.2_200.dat pdb_0.2_200.dat`

## Use of save_png.cpp

save_png.cpp skips the gtk plotter window and directly prints a png image showing the yz view of the beam through the scanner. The syntax to run save_png.cpp is 

`./save_pgn.run 0.2 200`

Remember to set the y and z range for the png to print in the cpp file.

## Automatic scan sequence using emittance_scan.sh

Use the shell script to run through scanning the y position and deflection voltage of the emittance scanner. The scans are in range `(seq start step stop)` for both parameters defined in the shell script.

The shell script supports continuing scan after breaks. The condition to determine which scan steps to skip is based on the existence of png files. If filenames match, the point will be skipped. It is optional to run the shell script with continue or restart.
