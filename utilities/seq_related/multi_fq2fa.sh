#!/bin/bash

for fq in $@; do
  seqtk seq -A $fq
done
