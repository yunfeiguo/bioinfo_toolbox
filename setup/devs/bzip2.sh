wget -O- http://www.bzip.org/1.0.6/bzip2-1.0.6.tar.gz | tar zxvf -
pushd bzip2-1.0.6
make
make install PREFIX=$HOME/Downloads/bzip2-1.0.6_install
popd
