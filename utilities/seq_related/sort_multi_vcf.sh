#!/bin/bash
echo "sorting and merging vcfs...
header from $1 will be used" >& 2
cat $1 | grep "^#"; cat $@ | grep -v "^#" | sort -k1,1d -k2,2n
