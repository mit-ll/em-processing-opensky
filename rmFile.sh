#!/bin/bash
# Copyright 2018 - 2020, MIT Lincoln Laboratory
# SPDX-License-Identifier: BSD-2-Clause

# $1 : Name of file to be deleted
#
# This script can be called directly or by a run script

# File to be removed
inFile=$1

# Display status to screen
echo "Removing: $inFile"

# Remove
rm -f $inFile
