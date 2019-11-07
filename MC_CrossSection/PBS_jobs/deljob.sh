#!/bin/sh

for id in {9612179..9612630}
do

echo qdel $id
qdel $id
done
