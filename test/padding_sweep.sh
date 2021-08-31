#!/bin/sh

for ptype in "global" "block" "edge"
do
    for filltype in "avg" "min" "max"
    do
        sh ./sweep_eb.sh $ptype $filltype
    done
done
