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
curl -L https://github.com/samtools/samtools/releases/download/1.4.1/samtools-1.4.1.tar.bz2 | tar jxvf -
pushd samtools-1.4.1
./configure --prefix $PWD
make -j
popd

#bedtools
wget -O - https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz | tar zxv
pushd bedtools2/
make -j
ln -sf $PWD/bin/* $HOME/bin
popd
#seqtk
git clone git@github.com:lh3/seqtk.git
cd seqtk; make
ln -s $PWD/seqtk $HOME/bin
cd ..
#freebayes
#htslib (contains tabix, bgzip)
git clone https://github.com/samtools/htslib.git
pushd htslib/
./configure --disable-lzma --disable-bz2
make -j
ln -sf $PWD/tabix ~/bin/
ln -sf $PWD/bgzip ~/bin
popd


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


#UCSC utilities
mkdir -p $HOME/Downloads/UCSC_util/
pushd $HOME/Downloads/UCSC_util/
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod +x *

mkdir -p chains
pushd chains
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
gunzip *.gz
popd

ln -sf $HOME/Downloads/UCSC_util/* $HOME/bin
popd


#blast
wget -O- ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.6.0+-x64-linux.tar.gz | tar zxvf -
ln -sf $PWD/ncbi-blast-2.6.0+/bin/*blast* $HOME/bin
