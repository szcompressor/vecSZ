prog='../bin/vecsz'
datadir='../../../datasets'

ptype=$1
filltype=$2

for bs in 8 16 32 64
do
for eb in 0.1 0.01 0.001 0.0001 0.00001 0.000001 0.0000001
do
    echo "PARAMS|Application=CESM|Dataset=CLDHGH_1_1800_3600|EB=$eb|Block=$bs|Padding=$ptype,$filltype"
	$prog -r -V -t f32 -D cesm -i $datadir/cesm/CLDHGH_1_1800_3600.f32 -e $eb -b $bs -p $ptype,$filltype
    echo "PARAMS|Application=Hurricane|Dataset=CLOUDf48|EB=$eb|Block=$bs|Padding=$ptype,$filltype"
	$prog -r -V -t f32 -D hurricane -i $datadir/hurr/CLOUDf48.bin.f32 -e $eb -b $bs -p $ptype,$filltype
done
for eb in 10 1 0.1 0.001
do
    echo "PARAMS|Application=NYX|Dataset=dark_matter_density|EB=$eb|Block=$bs|Padding=$ptype,$filltype"
	$prog -r -V -t f32 -D nyx -i $datadir/nyx/dark_matter_density.dat -e $eb -b $bs -p $ptype,$filltype
done
done
