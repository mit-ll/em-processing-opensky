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

if [ -z "$2" ]; then
    suffix="csv"
else
    suffix=$2
fi

outName=$AEM_DIR_OPENSKY/output/files_$suffix$(echo "$inDir" | sed 's/\//_/g').txt

echo "Output file: $outName"

find "$(cd $inDir ; pwd)" -name *$suffix > $outName

### LLSC
nproc=512
mapper="rmFile.sh"
dirListName=$outName
currDir=`pwd`

echo "LLMapReduce --input=$dirListName --output=$currDir
--mapper=$mapper --np=$nproc --keep true"

## Submit a job to the Grid to remove the directories in the file
LLMapReduce --input=$dirListName \
    --output=$currDir \
    --mapper=$mapper \
    --np=$nproc \
    --keep true
