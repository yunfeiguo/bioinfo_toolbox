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

#when need to run or compile programs dependent on gcc-5.4, just load this file by `source load_gcc-5.4.sh`
echo 'GCC_5_4_ROOT=$HOME/Downloads/gcc-5.4.0_install
export PATH=$GCC_5_4_ROOT/bin:$PATH
export CPPFLAGS="-I$GCC_5_4_ROOT/include $CPPFLAGS"
export LDFLAGS="-L$GCC_5_4_ROOT/lib $LDFLAGS"
export LD_LIBRARY_PATH=$GCC_5_4_ROOT/lib:$LD_LIBRARY_PATH' \
> $HOME/env/load_gcc-5.4.0.sh
