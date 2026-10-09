#!/bin/bash
make clean && make simu plot && ./simu.run && ./plot.run geom.dat epot.dat pdb.dat
