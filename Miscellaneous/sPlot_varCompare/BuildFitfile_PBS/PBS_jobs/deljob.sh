#!/bin/sh

for id in {6603613..6603645}
do

echo qdel $id
qdel $id
done
