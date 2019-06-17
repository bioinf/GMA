#!/bin/bash

genome="$(pwd)/""$1"
length=100
error=0.00
insertsize=0

function white_bold {
  echo -e "\033[1;37m$1\033[0m"
}
function white {
  echo -e "\033[37m$1\033[0m"
}
function stop {
  echo -e "\033[31mError:\n $1\033[0m\n";
  echo -e "\033[32mUsage:\n ./wrapper.sh genome.fasta -l 100 -s 0.00 -o 100\033[0m\n";
  exit 1;
}

while [ -n "$1" ]; do
  case "$1" in
    -l) length="$2"; shift;;
    -s) error="$2"; shift;;
    -o) insertsize="$2"; shift;;
    esac
  shift
done

if [[ $genome == "" ]]; then
	stop "Genome file not set"
fi
if [ ! -f $genome ]; then
  stop "Genome file not found – $genome"
fi

white "------------------------------------------------------------------------"
white_bold "GMA Wrapper"
white "  Genome: $genome"
white "  Length: $length"
white "  Error: $error"
white "  Insert Size: $insertsize"
white "------------------------------------------------------------------------"

RESULT=$(pwd)"/results_"$(date +'%m.%d.%Y_%H.%M.%S')
GMA_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
WORK_DIR="$(dirname $genome)"

mkdir $RESULT

"$GMA_DIR"/recordmem.sh "$RESULT/memory_used.log" &
recorder_pid=$!

trap 'kill -9 $recorder_pid' SIGINT

# Make INDEX
"$GMA_DIR"/bin/bwa index -a bwtsw "$genome"

# Preproduction
cd "$WORK_DIR"
"$GMA_DIR"/bin/prepro.chr.py "$genome"

# Production
PPD="$genome"".ppd"

cd $RESULT

i=0; d=0
cat $PPD | "$GMA_DIR"/bin/mapper runall -l $length -s $error -o $insertsize -q A -i $i -d $d -t 20 -f ref.fa -b 70 -x "$genome" -p $GMA_DIR/bin 1> map.txt

white_bold "GMA: Reducer runned"
cat map.txt | sort > mapsort.txt
cat mapsort.txt | "$GMA_DIR"/bin/reducer analyzer -l $length -t 20 -o $insertsize 1> mapred.txt 2> log

white_bold "Done"

kill $recorder_pid
