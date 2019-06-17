#!/bin/bash

logfile="$1"

function used {
  sec=`date +%s`
  echo -e "$sec\t"$(ps axo user:20,rss | grep $USER | awk '{ sum += $2 } END { printf "%.2f Mb\n", sum/1024; }')
}

while [ "TRUE" ]; do
  used >> $logfile
  sleep 1
done
