#!/bin/sh

set -eu

while read line; do
  bfn=`basename $line .o`
  \ls ${bfn}.F*
done

exit 0

