#!/bin/bash

#source setenv_papi.sh

echo "Running CESM Tests..."
export OMP_NUM_THREADS=32
for i in {1..10}
do
	for vers in pszO3_8 pszO3_16 pszO3_32 pszO3_64 pszvec_8 pszvec_16 pszvec_32 pszvec_64 pfetch_8 pfetch_16 pfetch_32 pfetch_64 pfetch2_8 pfetch2_16 pfetch2_32 pfetch2_64 pfetch4_8 pfetch4_16 pfetch4_32 pfetch4_64 pfetch8_8 pfetch8_16 pfetch8_32 pfetch8_64 pfetch32_8 pfetch32_16 pfetch32_32 pfetch32_64 pfetch32_8 pfetch32_16 pfetch32_32 pfetch32_64 pszomp_8 pszomp_16 pszomp_32 pszomp_64
	do
		for dset in CLDHGH_1_1800_3600.f32 CLDMED_1_1800_3600.f32 CLDLOW_1_1800_3600.f32
		do
			grep -q "pszomp" pomp || taskset 0x1 ./build/2D/$vers abs 1 -5 yesblk dq cesm /mydata/cesm/$dset
			grep -q "pszomp" pomp && ./build/2D/$vers abs 1 -5 yesblk dq cesm /mydata/cesm/$dset
		done
	done
done

echo "Running EXAFEL Tests..."
export OMP_NUM_THREADS=32
for i in {1..10}
do
	for vers in pszO3_8 pszO3_16 pszO3_32 pszO3_64 pszvec_8 pszvec_16 pszvec_32 pszvec_64 pfetch_8 pfetch_16 pfetch_32 pfetch_64 pfetch2_8 pfetch2_16 pfetch2_32 pfetch2_64 pfetch4_8 pfetch4_16 pfetch4_32 pfetch4_64 pfetch8_8 pfetch8_16 pfetch8_32 pfetch8_64 pfetch32_8 pfetch32_16 pfetch32_32 pfetch32_64 pfetch32_8 pfetch32_16 pfetch32_32 pfetch32_64 pszomp_8 pszomp_16 pszomp_32 pszomp_64	
	do

		for dset in smd-cxif5315-r129-dark.u16  smd-cxif5315-r169-calib-fde.u16  smd-cxif5315-r169-raw.u16
		do
			grep -q "pszomp" pomp || taskset 0x1 ./build/2D/$vers abs 1 -5 yesblk dq exafeldemo /mydata/exafel/$dset
			grep -q "pszomp" pomp && ./build/2D/$vers abs 1 -5 yesblk dq exafeldemo /mydata/exafel/$dset
		done
	done
done
