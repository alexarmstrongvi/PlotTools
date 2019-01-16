#!/bin/bash
#echo "TESTING :: BASH_SOURCE = ${BASH_SOURCE}"
if [ $0 != "-bash" ]; then
    echo "ERROR :: Must source directory"
else
    prog_dir="$(cd "$(dirname "$(dirname "$BASH_SOURCE")")"; pwd -P)"
    export PYTHONPATH="$PYTHONPATH:${prog_dir}"
    echo "$prog_dir added to PYTHONPATH"
    echo ">> echo PYTHONPATH"
    echo "."
    echo "."
    echo $PYTHONPATH | tr ":" "\n" | tail -n 3
fi
