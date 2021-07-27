prog='../bin/vecsz'
datadir='../../../datasets'

for bs in 8 16 32 64
do
for eb in 0.1 0.01 0.001 0.0001 0.00001 0.000001 0.0000001
do
	$prog -z -V -t f32 -D hacc -i $datadir/hacc/xx.f32 -e $eb -b $bs
	$prog -z -V -t f32 -D cesm -i $datadir/cesm/CLDHGH_1_1800_3600.f32 -e $eb -b $bs
	$prog -z -V -t f32 -D nyx -i $datadir/nyx/temperature.dat -e $eb -b $bs
	$prog -z -V -t f32 -D hurricane -i $datadir/hurr/CLOUDf48.bin.f32 -e $eb -b $bs
done
done
