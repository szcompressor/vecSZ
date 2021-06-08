#!/bin/bash

echo "Running HACC Tests..."
for dset in vx.f32 vy.f32 vz.f32 xx.f32 yy.f32 zz.f32
do
	for i in {1..10}
	do
		taskset 0x1 sz -z -f -M ABS -A 1E-4 -i /mydata/hacc/$dset -1 280953867
	done
done

echo "Running CESM Tests..."
for dset in CLDHGH_1_1800_3600.f32 CLDMED_1_1800_3600.f32 CLDLOW_1_1800_3600.f32
do
	for i in {1..10}
	do
		taskset 0x1 sz -z -f -M ABS -A 1E-5 -i /mydata/cesm/$dset -2 3600 1800
	done
done

echo "Running NYX Tests..."
for dset in baryon_density.dat dark_matter_density.dat temperature.dat velocity_x.dat velocity_y.dat velocity_z.dat
do
	for i in {1..10}
	do
		taskset 0x1 sz -z -f -M ABS -A 1E-4 -i /mydata/nyx/$dset -3 512 512 512
	done
done

echo "Running QMC Tests..."
for i in {1..10}
do
	taskset 0x1 sz -z -f -M ABS -A 1E-4 -i /mydata/qmc/einspline_288_115_69_69.pre.f32 -3 69 69 33120
done
for i in {1..10}
do
	taskset 0x1 sz -z -f -M ABS -A 1E-4 -i /mydata/qmc/einspline_115_69_69_288.f32 -3 288 69 7935
done

echo "Running Hurricane Tests..."
for dset in CLOUDf48.bin.f32 PRECIPf48.bin.f32 Pf48.bin.f32 QCLOUDf48.bin.f32 QGRAUPf48.bin.f32 QICEf48.bin.f32 QRAINf48.bin.f32 QSNOWf48.bin.f32 QVAPORf48.bin.f32 TCf48.bin.f32 Uf48.bin.f32 Vf48.bin.f32 Wf48.bin.f32
do
	for i in {1..10}
	do
		taskset 0x1 sz -z -f -M ABS -A 1E-4 -i /mydata/hurr/$dset -3 500 500 100
	done
done
