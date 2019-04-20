wget -O- https://download.sourceforge.net/libpng/libpng-1.6.32.tar.gz | tar zxvf -
pushd libpng-1.6.32
./configure --prefix=$HOME/Downloads/libpng-1.6.32_install
make
make install
