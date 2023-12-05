#!/bin/bash

# Check if the number of arguments is not between 1 and 3
if [ "$#" -lt 1 ] || [ "$#" -gt 3 ]; then
    echo "You must enter between 1 and 3 arguments: 'ordered', 'static', and/or 'wave'"
    exit 1
fi

# Loop over each argument
for dir in "$@"
do
    # Check if the argument is one of the allowed ones
    if [ "$dir" = "ordered" ] || [ "$dir" = "static" ] || [ "$dir" = "wave" ]; then
        # Delete all files in the directory
        rm -rf ./$dir/*
    else
        echo "Invalid argument: $dir. Allowed arguments are 'ordered', 'static', and 'wave'"
        exit 1
    fi
done

# run it with: 
# bash delete_files.sh ordered static wave
# to delete all files in the ordered, static, and wave directories