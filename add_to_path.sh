#!/bin/bash
#echo "TESTING :: BASH_SOURCE = ${BASH_SOURCE}"
if [ $0 != "-bash" ]; then
    echo "ERROR :: Must source directory"
else
    full_path="$(cd "$(dirname "$BASH_SOURCE")" && pwd -P)/$(basename "$1")"
    prog_dir="$(cd "$(dirname "$full_path")"; pwd -P)"
    export PYTHONPATH="$PYTHONPATH:${prog_dir}"
    #echo "$prog_dir added to PYTHONPATH"
    #echo 
    #echo "See for yourself..."
    #echo "Running: echo \$PYTHONPATH | tr \":\" \"\\n\" | tail -n 3"
    #echo $PYTHONPATH | tr ":" "\n" | tail -n 3
fi
