#!/bin/bash

MEM_REGEX="(0.125m|0.25m|0.375m|0.5m|0.625m|0.75m|0.875m|1m)"
SKETCH_REGEX="(ES|FR|NS|NZE|PR|UM|BS_ES|BS_FR|BS_NS|BS_NZE|BS_PR|BS_UM)"

# Helper function
DisplayHelp()
{
   echo "OVERVIEW: Run per-flow measurement"
   echo
   echo "USAGE: $0 -m mem-size -s sketch-name [-v|-h]"
   echo
   echo "OPTIONS:"
   echo "  m     Specify the memory budget $MEM_REGEX"
   echo "  s     Specify the sketch framework $SKETCH_REGEX"
   echo "  v     Verbose mode"
   echo "  h     Print this help message"
   echo
}

# Error
function Die() {
    echo >&2 "$@"
    exit 1
}

# Arguments
verbose=0
memory=""
sketch=""
## Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -m|--memory)
          memory="$2";
          shift
        ;;
        -s|--sketch)
          sketch="$2"; 
          shift
        ;;
        -v|--verbose)
          verbose=1
        ;;
        -h|--help)
          DisplayHelp ;
          exit 0
        ;;
        *)
          Die "ERROR: Unknown parameter passed: $1";
        ;;
    esac
    shift
done
## Check arguments' validity
echo "$memory" | grep -E -q "^$MEM_REGEX" \
  || Die "ERROR: Invalid memory size, which should be in $MEM_REGEX but got \"$memory\" instead"
echo "$sketch" | grep -E -q "^$SKETCH_REGEX" \
  || Die "ERROR: Invalid sketch name, which should be in $SKETCH_REGEX but got \"$sketch\" instead"


# Run sketch framework
## Construct path to the config file
function PathToConfig () {
  if ( echo "$1" | grep -E -q '^BS.*' ); then
    prefix="${sketch:3}"
    echo "../config/PerFlowMeas/Optimized/$prefix/$memory.toml"
  else
    prefix="$sketch"
    echo "../config/PerFlowMeas/Raw/$prefix/$memory.toml"
  fi
}
## Verbose mode
function ConfigMessage () {
  echo "Reading sketch config from \"$1\":"
  echo
  cat $1
  echo
  echo "======"
}
## Parse Output
function ParseOutput () {
  if (echo "$1" | grep -E -q '.*FR$'); then # FR / BS_FR
    arg="${2#*DECODE COST}" # Remove some head logs
    MEM=$(echo "$arg" | awk -F "Footprint:" '{print $2}' | tr -d [:space:])
    MEM=$(expr "$MEM" : '\([^B]*B\)') # Remove multiple displays
    ARE=$(echo "$arg" | awk -F "Decode ARE:" '{print $2}' | tr -d [:space:])
    local DecodeRatio=$(echo "$arg" | awk -F "Decode Ratio:" '{print $2}' | tr -d [:space:%])
    local DecodePODF=$(echo "$arg" | awk -F "Decode <=0.1%:" '{print $2}' | tr -d [:space:%])
    FR=$(python -c "print('{val:f}%'.format(val = $DecodeRatio * $DecodePODF / 100))")
  else # Others
    arg="${2#*DataSet}" # Remove some head logs
    MEM=$(echo "$arg" | awk -F "Footprint:" '{print $2}' | tr -d [:space:])
    MEM=$(expr "$MEM" : '\([^B]*B\)') # Remove multiple displays
    ARE=$(echo "$arg" | awk -F "Query ARE:" '{print $2}' | tr -d [:space:])
    FR=$(echo "$arg" | awk -F "<=0.1%:" '{print $2}' | tr -d [:space:])
  fi
  echo "Mem: $MEM, ARE: $ARE, FR: $FR"
}
## Check the executable and config file exist
path=$(PathToConfig $sketch)
if [ ! -f $path ]; then
  Die "ERROR: Config file \"$path\" does not exist"
fi
if [ ! -f $sketch ]; then
  Die "ERROR: Sketch framework \"$sketch\" does not exist"
fi
if [ ! -x $sketch ]; then
  Die "ERROR: Sketch framework \"$sketch\" is not executable"
fi
if [ $verbose -gt 0 ]; then
  ConfigMessage $path
fi
## Run sketch
output=$(./$sketch -c $path 2>/dev/null)
if [ $verbose -gt 0 ]; then
  echo "$output"
else
  ParseOutput "$sketch" "$output"
fi
