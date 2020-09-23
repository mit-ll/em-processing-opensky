#!/bin/bash
#
# $1 : Name of directory to be deleted
#
# This script can be called directly or by a run script

# zip man 
# -j --junk-paths Store just the name of a saved file (junk the path), and do not store directory names. By default, zip will store the full path (relative to the current directory).
# -g --grow Grow (append to) the specified zip archive, instead of creating a new one. If this operation fails, zip attempts to restore the archive to its original state. If the restoration fails, the archive might become corrupted. This option is ignored when there's no existing archive or when at least one archive member must be updated or deleted.
# -u --update Update existing entries if newer on the file system and add new files. If the archive does not exist issue warning then create a new archive.

# Directory to be archived
inDir=$1

# Archive to create
outFile=$(echo "$inDir" | sed "s/1_organize/2_archive/")
outFile=$outFile.zip

# Only do something if there are files
# https://unix.stackexchange.com/a/202352/1408
if [ ! -z "$(ls -A -- "$inDir")" ]; then
    # Archive
    if [ -f "$outFile" ]; then
        echo "Updating: $outFile"
        zip -j -r -u $outFile $inDir
    else
        echo "Archive does not exist, creating: $outFile"
        zip -j -r $outFile $inDir
    fi
else
    echo "Directory empty: $inDir"
fi
