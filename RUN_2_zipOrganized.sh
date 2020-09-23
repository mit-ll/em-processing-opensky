#!/bin/bash
# Copyright 2018 - 2020, MIT Lincoln Laboratory
# SPDX-License-Identifier: BSD-2-Clause

# cleanupSubmitted.sh
#
# $1 : name of file containing list of directories to delete

# Directory for output of oranizeraw_1()
if [ -z "$1" ]; then
    inDir=$AEM_DIR_OPENSKY/output/1_organize
else
    inDir=$1
fi

# Find directory names
# https://askubuntu.com/a/444554/244714
find "$(cd $inDir ; pwd)" -mindepth 1 -maxdepth 1 -type d> $AEM_DIR_OPENSKY/output/2_dirArchiveDepth1.txt
find "$(cd $inDir ; pwd)" -mindepth 2 -maxdepth 2 -type d> $AEM_DIR_OPENSKY/output/2_dirArchiveDepth2.txt
find "$(cd $inDir ; pwd)" -mindepth 3 -maxdepth 3 -type d> $AEM_DIR_OPENSKY/output/2_dirArchiveDepth3.txt
find "$(cd $inDir ; pwd)" -mindepth 4 -maxdepth 4 -type d> $AEM_DIR_OPENSKY/output/2_dirArchiveDepth4.txt

# Replicate 1_organize directories with depth <= 3 in 2_archive directory
# https://stackoverflow.com/a/1521498/363829
# https://stackoverflow.com/a/22957485/363829
while read d; do
    newDir=$(echo "$d" | sed "s/1_organize/2_archive/")
    if [ ! -d "$newDir" ]; then
        echo $newDir
        mkdir -p $newDir
    fi
done < $AEM_DIR_OPENSKY/output/2_dirArchiveDepth3.txt

### Serial - Depth 3
#while read d; do
#    newDir="${d//1_organize/2_archive}"
#    echo $newDir
#    mkdir -p $newDir
#done < $AEM_DIR_OPENSKY/output/2_dirArchiveDepth3.txt

### LLSC - Depth 3
nproc=512
mapper="zipOrganized_2.sh"
dirListName=$AEM_DIR_OPENSKY/output/2_dirArchiveDepth4.txt
currDir=`pwd`

echo "LLMapReduce --input=$dirListName --output=$currDir
--mapper=$mapper --np=$nproc --keep true"

## Submit a job to the Grid to remove the directories in the file
LLMapReduce --input=$dirListName \
    --output=$currDir \
    --mapper=$mapper \
    --np=$nproc \
    --keep true
