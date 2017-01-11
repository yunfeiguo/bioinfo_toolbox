#!/bin/bash
set -ex
#bwa
git clone https://github.com/lh3/bwa.git
cd bwa; make
ln -sf $PWD/bwa $HOME/bin/
#bamtools
git clone git@github.com:pezmaster31/bamtools.git
cd bamtools
mkdir build
cd build
cmake ..
make
cd ..
ln -sf $PWD/bin/bamtools ~/bin/
#samtools
#bedtools
#seqtk
git clone https://github.com/lh3/seqtk.git;
cd seqtk; make
ln -s $PWD/seqtk $HOME/bin
cd ..
#freebayes
#tabix
git clone git@github.com:samtools/tabix.git
cd tabix/
make
ln -s $PWD/tabix ~/bin/
ln -s $PWD/bgzip ~/bin/
