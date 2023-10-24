#!/bin/bash

SKETCH_REGEX="(ES|FR|NS|NZE|PR|UM|BS_ES|BS_FR|BS_NS|BS_NZE|BS_PR|BS_UM)"

# Helper function
DisplayHelp()
{
   echo "OVERVIEW: Helper of per-flow measurement"
   echo "          Run sketch on all memory settings"
   echo
   echo "USAGE: $0 -o output-file -s sketch-name [-p proc-number] [-h]"
   echo
   echo "OPTIONS:"
   echo "  s     Specify the sketch framework $SKETCH_REGEX"
   echo "  p     Specify the number of processes to run in parallel (default = 1)"
   echo "  o     Specify the output file"
   echo "  h     Print this help message"
   echo
}

# Error
function Die() {
    echo >&2 "$@"
    exit 1
}

# Arguments
process=1
sketch=""
file=""
## Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -s|--sketch)
          sketch="$2"; 
          shift
        ;;
        -p|--process)
          process="$2";
          shift
        ;;
        -o|--output)
          file="$2";
          shift
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
echo "$sketch" | grep -E -q "^$SKETCH_REGEX" \
  || Die "ERROR: Invalid sketch name, which should be in $SKETCH_REGEX but got \"$sketch\" instead"
echo "$process" | grep -E -q '^[1-9][0-9]*$' \
  || Die "ERROR: Invalid number of processes, got \"$process\""
if [ -z "$file" ]; then
  Die "ERROR: No output file specified (Please use -o)"
fi
if ( ! touch "$file" ); then
  Die "ERROR: Touch file \"$file\" failed"
fi

# Run single sketch
function RunSingleSketch () {
  echo "$1" > "$file"
  for memory in "0.125m" "0.25m" "0.375m" "0.5m" "0.625m" "0.75m" "0.875m" "1m"; do 
    ((i=i%"$2")); ((i++==0)) && wait
    # Warning: Potential contention
    #   Fortunately each process appends a single line to the shared file,
    # and there are at most 8 processes. This implementation is for simp-
    # licity only.
    ./PerFlowMeas.sh -s "$1" -m $memory >> "$file" &
  done
  wait
}

RunSingleSketch $sketch $process