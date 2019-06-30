if [ "$BASH_SOURCE" == "$(basename "$BASH_SOURCE")" ]; then
    # sourcing file in atlasrootstyle

    full_path="$(dirname "$(pwd -P)")"
else
    # sourcing file from outside atlasrootstyle
    gparent_dir="$(dirname "$(dirname "$BASH_SOURCE")")"
    full_path="$(cd $gparent_dir && pwd -P)"
fi
export PYTHONPATH="$PYTHONPATH:${full_path}"
echo "$full_path added to PYTHONPATH. Now you can run \"from PlotTools import ...\""
