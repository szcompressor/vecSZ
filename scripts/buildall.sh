#!/bin/bash
[ -f "./build" ] || mkdir -p build/{1D,2D,3D}

make clean
for dim in 1D 2D 3D
do
	sed -i "s/DIM= -D\_.*/DIM= -D\_$dim/" Makefile
	sed -i "s/DIMDIR?=.*/DIMDIR?=$dim/" Makefile
	make -s vec
	make -s O3
	make -s pfetch
	make -s pfetch2
	make -s pfetch4
	make -s pfetch8
	make -s pfetch16
	make -s ompsz
	make -s auto
done
