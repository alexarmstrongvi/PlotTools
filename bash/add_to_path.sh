#!/bin/bash
if [ $0 != "-bash" ]; then
    echo "ERROR :: Must source file"
else
    full_path="$(cd "$(dirname "$(dirname "$BASH_SOURCE")")" && pwd -P)"
    prog_dir="$(cd "$(dirname "$full_path")"; pwd -P)"
    export PYTHONPATH="$PYTHONPATH:${prog_dir}"
    echo "${prog_dir} added to PYTHONPATH"
fi
