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
git clone git@github.com:lh3/seqtk.git
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



#ngmlr, make sure to have at least gcc 5.3
git clone https://github.com/philres/ngmlr.git
pushd ngmlr
mkdir build
pushd build
cmake ..
#modify LIBRARY_PATH as necessary because the compiler looks for static libraries (lib*.a) during compilation
#changing LD_LIBRARY_PATH does not help as it is used for shared libraries (after compilation is done)
make -j
#cmake does not seem to have good support for make -j
make
ln -s $PWD/
popd
popd
