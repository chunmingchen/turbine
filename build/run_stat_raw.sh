#!/bin/bash
dir=/data/turbine_Stg/zDIR.P3D.rel.6201-11001/
for f in $dir/*.list
do
    id=$(echo $f | sed -e 's/.*\///g')
    id=$(echo $id | sed -e 's/\.list//g')
    cmd="./stat_raw $f $id.vtp"
    echo $cmd
    $cmd
done

