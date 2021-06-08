#!/bin/bash

#source setenv_papi.sh
echo "Running HACC Tests..."
export OMP_NUM_THREADS=32
for i in {1..5}
do
	for nt in 1 2 4 8 16 32 64
	do
	for vers in pszomp_8 pszomp_16 pszomp_32 pszomp_64
	do
		for dset in vx.f32 vy.f32 vz.f32 xx.f32 yy.f32 zz.f32
		do
			export OMP_NUM_THREADS=$nt
			export OMP_PLACES=cores
			export OMP_PROC_BIND=close
			echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
			 ./build/1D/$vers abs 1 -4 yesblk dq hacc /mydata/hacc/$dset
		done
	done
	done
done

echo "Running CESM Tests..."
for i in {1..5}
do
	for nt in 1 2 4 8 16 32 64
	do
	for vers in pszomp_8 pszomp_16 pszomp_32 pszomp_64
	do for dset in CLDHGH_1_1800_3600.f32 CLDMED_1_1800_3600.f32 CLDLOW_1_1800_3600.f32
		do
			export OMP_NUM_THREADS=$nt
			export OMP_PLACES=cores
			export OMP_PROC_BIND=close
			echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
			 ./build/2D/$vers abs 1 -5 yesblk dq cesm /mydata/cesm/$dset
		done
	done
	done
done

echo "Running NYX Tests..."
for i in {1..5}
do
	for nt in 1 2 4 8 16 32 64
	do
	for vers in pszomp_8 pszomp_16 pszomp_32 pszomp_64
	do
		for dset in baryon_density.dat dark_matter_density.dat temperature.dat velocity_x.dat velocity_y.dat velocity_z.dat
		do
			export OMP_NUM_THREADS=$nt
			export OMP_PLACES=cores
			export OMP_PROC_BIND=close
			echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
			 ./build/3D/$vers abs 1 -4 yesblk dq nyx /mydata/nyx/$dset
		done
	done
	done
done

echo "Running QMC Tests..."
for i in {1..5}
do
	for nt in 1 2 4 8 16 32 64
	do
	for vers in pszomp_8 pszomp_16 pszomp_32 pszomp_64
	do
		for dset in einspline_288_115_69_69.pre.f32 einspline_115_69_69_288.f32
		do
			export OMP_NUM_THREADS=$nt
			export OMP_PLACES=cores
			export OMP_PROC_BIND=close
			echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
			 ./build/3D/$vers abs 1 -4 yesblk dq qmcpre /mydata/qmc/$dset
		done
	done
	done
done

echo "Running Hurricane Tests..."
for i in {1..5}
do
	for nt in 1 2 4 8 16 32 64
	do
	for vers in pszomp_8 pszomp_16 pszomp_32 pszomp_64
	do
		for dset in CLOUDf48.bin.f32 PRECIPf48.bin.f32 Pf48.bin.f32 QCLOUDf48.bin.f32 QGRAUPf48.bin.f32 QICEf48.bin.f32 QRAINf48.bin.f32 QSNOWf48.bin.f32 QVAPORf48.bin.f32 TCf48.bin.f32 Uf48.bin.f32 Vf48.bin.f32 Wf48.bin.f32
		do
			export OMP_NUM_THREADS=$nt
			export OMP_PLACES=cores
			export OMP_PROC_BIND=close
			echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
			 ./build/3D/$vers abs 1 -4 yesblk dq hurricane /mydata/hurr/$dset
		done
	done
	done
done
