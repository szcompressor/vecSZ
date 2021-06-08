#!/bin/bash

#source setenv_papi.sh

echo "Running NYX Tests..."
export OMP_NUM_THREADS=32
for i in {1..10}
do
	for vers in pszO3_8 pszO3_16 pszO3_32 pszO3_64 pszvec_8 pszvec_16 pszvec_32 pszvec_64 pfetch_8 pfetch_16 pfetch_32 pfetch_64 pfetch2_8 pfetch2_16 pfetch2_32 pfetch2_64 pfetch4_8 pfetch4_16 pfetch4_32 pfetch4_64 pfetch8_8 pfetch8_16 pfetch8_32 pfetch8_64 pfetch32_8 pfetch32_16 pfetch32_32 pfetch32_64 pfetch32_8 pfetch32_16 pfetch32_32 pfetch32_64 pszomp_8 pszomp_16 pszomp_32 pszomp_64
	do
		for dset in baryon_density.dat dark_matter_density.dat temperature.dat velocity_x.dat velocity_y.dat velocity_z.dat
		do
			grep -q "pszomp" pomp || taskset 0x1 ./build/3D/$vers abs 1 -4 yesblk dq nyx /mydata/nyx/$dset
			grep -q "pszomp" pomp && ./build/3D/$vers abs 1 -4 yesblk dq nyx /mydata/nyx/$dset
		done
	done
done

echo "Running QMC Tests..."
export OMP_NUM_THREADS=32
for i in {1..10}
do
	for vers in pszO3_8 pszO3_16 pszO3_32 pszO3_64 pszvec_8 pszvec_16 pszvec_32 pszvec_64 pfetch_8 pfetch_16 pfetch_32 pfetch_64 pfetch2_8 pfetch2_16 pfetch2_32 pfetch2_64 pfetch4_8 pfetch4_16 pfetch4_32 pfetch4_64 pfetch8_8 pfetch8_16 pfetch8_32 pfetch8_64 pfetch32_8 pfetch32_16 pfetch32_32 pfetch32_64 pfetch32_8 pfetch32_16 pfetch32_32 pfetch32_64 pszomp_8 pszomp_16 pszomp_32 pszomp_64
	do
		for dset in einspline_288_115_69_69.pre.f32 
		do
			grep -q "pszomp" pomp || taskset 0x1 ./build/3D/$vers abs 1 -4 yesblk dq qmcpre /mydata/qmc/$dset
			grep -q "pszomp" pomp && ./build/3D/$vers abs 1 -4 yesblk dq qmcpre /mydata/qmc/$dset
		done
	done
done
