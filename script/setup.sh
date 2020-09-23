#!/bin/sh
# Copyright 2018 - 2020, MIT Lincoln Laboratory
# SPDX-License-Identifier: BSD-2-Clause

# script/setup:
# Set up application for the first time after cloning, 
# or set it back to the initial first unused state.

####### OpenSky Monday State Data

URL_OSN_STATES='https://opensky-network.org/datasets/states'

# Directory for raw data
if [ -z "$3" ]; then
    DIR_DOWNLOAD_ROOT=$AEM_DIR_OPENSKY/data
else
    DIR_DOWNLOAD_ROOT=$3
fi

# https://stackoverflow.com/a/28226339
# examples
#input_start=2020-03-16
#input_end=2020-06-01
input_start=$1
input_end=$2

# After this, startdate and enddate will be valid ISO 8601 dates,
# or the script will have aborted when it encountered unparseable data
# such as input_end=abcd
startdate=$(date -I -d "$input_start") || exit -1
enddate=$(date -I -d "$input_end")     || exit -1

d="$startdate"
while [ "$(date -d "$d" +%Y%m%d)" -lt "$(date -d "$enddate" +%Y%m%d)" ]; do
	echo $d
	DIR_DOWNLOAD=$DIR_DOWNLOAD_ROOT/$d
	mkdir -p $DIR_DOWNLOAD

	# Iterate over each hour
	# https://stackoverflow.com/a/8789815
	for hour in $(seq -f "%02g" 0 23)
	do
		URL_CURRENT=$URL_OSN_STATES/$d/$hour/states_$d-$hour.csv.tar
		FILE_DOWNLOAD=$DIR_DOWNLOAD/states_$d-$hour.csv.tar
		
		echo $FILE_DOWNLOAD

		# Download file 
		wget $URL_CURRENT -O $FILE_DOWNLOAD

		# Extract tar
		tar xvf $FILE_DOWNLOAD -C $DIR_DOWNLOAD

		# Extract again
		gunzip $DIR_DOWNLOAD/states_$d-$hour.csv.gz
	done

	# Advance day
	d=$(date -I -d "$d + 7 day")
done