#!/bin/bash

#source setenv_papi.sh
export OMP_NUM_THREADS=32
echo "Running NYX Tests..."
for i in {1..10}
do
	for nt in 1 2 4 8 16 32 64
	do
		for dset in baryon_density.dat dark_matter_density.dat temperature.dat velocity_x.dat velocity_y.dat velocity_z.dat
		do
			export OMP_NUM_THREADS=$nt
			echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
		 	sz -k -z -f -M ABS -A 1E-4 -i /mydata/nyx/$dset -3 512 512 512
		done
	done
done

echo "Running QMC Tests..."
for i in {1..10}
do
	for nt in 1 2 4 8 16 32 64
	do
			export OMP_NUM_THREADS=$nt
			echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
	 		sz -k -z -f -M ABS -A 1E-4 -i /mydata/qmc/einspline_288_115_69_69.pre.f32 -3 69 69 33120
	 		sz -k -z -f -M ABS -A 1E-4 -i /mydata/qmc/einspline_115_69_69_288.f32 -3 288 69 7935
	done
done

echo "Running Hurricane Tests..."
for i in {1..10}
do
	for nt in 1 2 4 8 16 32 64
	do
	for dset in CLOUDf48.bin.f32 PRECIPf48.bin.f32 Pf48.bin.f32 QCLOUDf48.bin.f32 QGRAUPf48.bin.f32 QICEf48.bin.f32 QRAINf48.bin.f32 QSNOWf48.bin.f32 QVAPORf48.bin.f32 TCf48.bin.f32 Uf48.bin.f32 Vf48.bin.f32 Wf48.bin.f32
	do
		export OMP_NUM_THREADS=$nt
		echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
		sz -k -z -f -M ABS -A 1E-4 -i /mydata/hurr/$dset -3 500 500 100
	done
	done
done
