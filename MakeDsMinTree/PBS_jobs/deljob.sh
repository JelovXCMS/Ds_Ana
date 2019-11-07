#!/bin/sh

for id in {13443515..13443929}
do

echo qdel $id
qdel $id
done
