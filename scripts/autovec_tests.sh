#!/bin/bash

cd /users/gadube/cuSZ/src/pSZ/

echo "Running CESM Tests..."
for niterations in 1 10 20 30  
do
for nsamples in 0.01 0.1 1 10 20 50
do
for dset in CLDHGH_1_1800_3600.f32 CLDMED_1_1800_3600.f32 CLDLOW_1_1800_3600.f32
do
	for i in {1..5}
	do
		taskset 0x1 ./build512/2D/pszvec_auto abs 1 -5 yesblk dq cesm /mydata/cesm/$dset $niterations $nsamples
	done
done
done
done

echo "Running NYX Tests..."
for niterations in 1 10 20 30  
do
for nsamples in 0.01 0.1 1 10 20 50
do
for dset in baryon_density.dat dark_matter_density.dat temperature.dat velocity_x.dat velocity_y.dat velocity_z.dat
do
	for i in {1..5}
	do
		taskset 0x1 ./build512/3D/pszvec_auto abs 1 -4 yesblk dq nyx /mydata/nyx/$dset $niterations $nsamples
	done
done
done
done

echo "Running QMC Tests..."
for niterations in 1 10 20 30  
do
for nsamples in 0.01 0.1 1 10 20 50
do
for i in {1..5}
do
	taskset 0x1 ./build512/3D/pszvec_auto abs 1 -4 yesblk dq qmcpre /mydata/qmc/einspline_288_115_69_69.pre.f32
done
for i in {1..5}
do
	taskset 0x1 ./build512/3D/pszvec_auto abs 1 -4 yesblk dq qmc /mydata/qmc/einspline_115_69_69_288.f32
done
done
done

echo "Running Hurricane Tests..."
for niterations in 1 10 20 30  
do
for nsamples in 0.01 0.1 1 10 20 50
do
for dset in CLOUDf48.bin.f32 PRECIPf48.bin.f32 Pf48.bin.f32 QCLOUDf48.bin.f32 QGRAUPf48.bin.f32 QICEf48.bin.f32 QRAINf48.bin.f32 QSNOWf48.bin.f32 QVAPORf48.bin.f32 TCf48.bin.f32 Uf48.bin.f32 Vf48.bin.f32 Wf48.bin.f32
do
	for i in {1..5}
	do
		taskset 0x1 ./build512/3D/pszvec_auto abs 1 -4 yesblk dq hurricane /mydata/hurr/$dset $niterations $nsamples
	done
done
done
done

echo "Running HACC Tests..."
for niterations in 1 10 20 30  
do
for nsamples in 0.01 0.1 1 10 20 50
do
for dset in vx.f32 vy.f32 vz.f32 xx.f32 yy.f32 zz.f32
do
	for i in {1..5}
	do
		taskset 0x1 ./build512/1D/pszvec_auto abs 1 -4 yesblk dq hacc /mydata/hacc/$dset $niterations $nsamples
	done
done
done
done
