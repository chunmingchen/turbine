#!/bin/bash
dir=/data/flow2/turbine_Stg/zDIR.P3D.rel.6201-11001/
for f in $dir/*.list.vti
do
    id=$(echo $f | sed -e 's/.*\///g')
    id=$(echo $id | sed -e 's/\.list.*//g')
    ./stat $f $id

done

