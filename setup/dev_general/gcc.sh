#download the tarball, unzip it, configure, build and install it
#create symbolic links in $HOME/bin
#write env loading scripts for future loading
wget -O- http://www.netgull.com/gcc/releases/gcc-5.4.0/gcc-5.4.0.tar.bz2 | tar jxvf -
pushd gcc-5.4.0/
./contrib/download_prerequisites
popd
mkdir -p gcc-5.4.0_build
pushd gcc-5.4.0_build
../gcc-5.4.0/configure --prefix=$HOME/Downloads/gcc-5.4.0_install \
  --enable-shared --enable-threads --with-system-zlib --enable-languages=c,c++,fortran,go --disable-multilib
  #make -j may have errors
  make && make install
ln -sf $HOME/Downloads/gcc-5.4.0_install/bin/gcc $HOME/bin
ln -sf $HOME/Downloads/gcc-5.4.0_install/bin/g++ $HOME/bin
popd

