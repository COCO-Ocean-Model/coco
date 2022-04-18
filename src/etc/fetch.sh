#!/bin/sh

set -eu

srcdir=$1
shift

for x in "$@"; do
    ln -sf "${srcdir}/$x" .
done
