mkdir -p $HOME/Downloads
mkdir -p $HOME/bin
cd $HOME/Downloads
NUM_THREADS=$( cat /proc/cpuinfo | grep -c '^processor')
USABLE_THREADS=$(expr $NUM_THREADS - 1)
#java
wget -O - --no-check-certificate --no-cookies --header "Cookie: oraclelicense=accept-securebackup-cookie" http://download.oracle.com/otn-pub/java/jdk/8u112-b15/jdk-8u112-linux-x64.tar.gz | tar zxvf -
ln -sf $PWD/jdk1.8.0_112/bin/* $HOME/bin

#glibc
wget https://ftp.gnu.org/gnu/libc/glibc-2.24.tar.gz
tar zxvf glibc-2.24.tar.gz 
mkdir glibc-2.24-objdir
pushd glibc-2.24-objdir/
../glibc-2.24/configure --prefix=$HOME/Downloads/glibc-2.24_install
make -j
make install
popd

#gcc5.4

enter the following into `load_gcc-5.4.sh`
####
#GCC_5_4_ROOT=$HOME/Downloads/gcc-5.4.0_install
#export PATH=$GCC_5_4_ROOT/bin:$PATH
#export CPPFLAGS="-I$GCC_5_4_ROOT/include $CPPFLAGS"
#export LDFLAGS="-L$GCC_5_4_ROOT/lib $LDFLAGS"
#export LD_LIBRARY_PATH=$GCC_5_4_ROOT/lib:$LD_LIBRARY_PATH
###
#when need to run or compile programs dependent on gcc-5.4, just load this file by `source load_gcc-5.4.sh`

#gcc, may a a few hours
wget http://www.netgull.com/gcc/releases/gcc-6.3.0/gcc-6.3.0.tar.gz
tar zxf gcc-6.3.0.tar.gz
pushd gcc-6.3.0/
./contrib/download_prerequisites
popd
mkdir gcc-6.3.0-objdir
pushd gcc-6.3.0-objdir/
../gcc-6.3.0/configure --prefix=$HOME/Downloads/gcc-6.3.0_install --enable-shared --enable-threads \
--with-system-zlib --enable-languages=c,c++,fortran,go --disable-multilib
#--with-gmp
#--with-mpfr
#--with-mpc
#--with-isl
make -j $USABLE_THREADS && make install
ln -sf $PWD/gcc $HOME/bin/
ln -sf $PWD/g++ $HOME/bin/
ln -sf $PWD/go $HOME/bin/
ln -sf $PWD/gfortran $HOME/bin/
popd
echo 'export LD_LIBRARY_PATH=$HOME/Downloads/gcc-6.3.0_install/lib:$LD_LIBRARY_PATH' >> $HOME/.bashrc
echo 'export LD_LIBRARY_PATH=$HOME/Downloads/gcc-6.3.0_install/lib64:$LD_LIBRARY_PATH' >> $HOME/.bashrc
echo 'export CC=$HOME/bin/gcc' >> $HOME/.bashrc
echo 'export CXX=$HOME/bin/g++' >> $HOME/.bashrc
#cmake
wget https://cmake.org/files/v3.7/cmake-3.7.2-Linux-x86_64.sh
bash cmake-3.7.2-Linux-x86_64.sh
ln -s $PWD/cmake-3.7.2-Linux-x86_64/bin/* $HOME/bin/


#zlib
wget http://www.zlib.net/zlib-1.2.11.tar.gz
tar zxvf zlib-1.2.11.tar.gz 
pushd  zlib-1.2.11/
./configure --prefix=$HOME/Downloads/zlib-1.2.11_install
make && make install
echo 'export LD_LIBRARY_PATH=$HOME/Downloads/zlib-1.2.11_install:$LD_LIBRARY_PATH' >> /home/guoy28/.bashrc
popd

#openssl (for https)
wget https://github.com/openssl/openssl/archive/OpenSSL_1_1_0e.tar.gz
tar zxvf OpenSSL_1_1_0e.tar.gz 
pushd openssl-OpenSSL_1_1_0e/
./config --prefix=$HOME/Downloads/openssl_1_1_0e_install && make -j && make install
popd

#libcurl
wget https://curl.haxx.se/download/curl-7.53.1.tar.gz
tar zxvf curl-7.53.1.tar.gz 
pushd curl-7.53.1/
./configure --prefix=$HOME/Downloads/curl-7.53.1_install --with-ssl=$HOME/Downloads/openssl_1_1_0e_install && make -j && make install
popd

#BLAS
wget http://github.com/xianyi/OpenBLAS/archive/v0.2.19.tar.gz
tar xzvf v0.2.19.tar.gz 
pushd OpenBLAS-0.2.19/
make -j
make PREFIX=$HOME/Downloads/OpenBLAS-0.2.19_install install
popd

#python +setuptools+pip

#curl -L https://prdownloads.sourceforge.net/tcl/tcl8.6.6-src.tar.gz | tar zxvf -
#pushd tcl8.6.6/unix/
#./configure --prefix=$HOME/Downloads/tcl8.6.6_install
#make -j && make install
#popd
#curl -L https://prdownloads.sourceforge.net/tcl/tk8.6.6-src.tar.gz | tar zxvf -
#pushd ../tk8.6.6/unix/
#./configure --prefix=$HOME/Downloads/tk8.6.6_install --with-tcl=$HOME/Downloads/tcl8.6.6/unix && make -j && make install
#popd

curl -L https://www.python.org/ftp/python/2.7.12/Python-2.7.12.tgz | tar zxvf -
pushd Python-2.7.12/
./configure  --prefix=$HOME/Downloads/python-2.7.12_install && make -j && make install
ln -sf $HOME/Downloads/python-2.7.12_install/bin/python $HOME/bin/
export PATH=$HOME/Downloads/python-2.7.12_install/bin:$PATH

wget https://github.com/pypa/setuptools/archive/v29.0.1.zip
unzip v29.0.1.zip 
pushd setuptools-29.0.1/
python bootstrap.py 
python setup.py install
popd

curl -L https://pypi.python.org/packages/11/b6/abcb525026a4be042b486df43905d6893fb04f05aac21c32c638e939e447/pip-9.0.1.tar.gz | tar zxvf - 
pushd pip-9.0.1/
python setup.py install
ln -sf $HOME/Downloads/python-2.7.12_install/bin/pip ~/bin/
popd

pip install numpy scipy matplotlib pysam pyvcf biopython

popd

#python3
wget -O - https://www.python.org/ftp/python/3.6.1/Python-3.6.1.tgz | tar zxf -
pushd Python-3.6.1/
./configure --prefix=$HOME/Downloads/Python-3.6.1_install
make -j && make install
ln -sf $HOME/Downloads/Python-3.6.1_install/bin/python3 $HOME/bin
ln -sf $HOME/Downloads/Python-3.6.1_install/bin/pip3 $HOME/bin

#perl
#boost
wget -O boost_1_63_0.tar.gz 'https://downloads.sourceforge.net/project/boost/boost/1.63.0/boost_1_63_0.tar.gz?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fboost%2Ffiles%2Fboost%2F1.63.0%2F&ts=1487717694&use_mirror=droneda
ta'
tar zxvf boost_1_63_0.tar.gz
pushd boost_1_63_0/
./bootstrap.sh 
./b2
popd

#xslt
if ! [ -x "$(command -v xsltproc)" ]; then
  wget ftp://xmlsoft.org/xslt/libxslt-1.1.29-rc1.tar.gz
  tar zxvf libxslt-1.1.29-rc1.tar.gz
  pushd libxslt-1.1.29/
  ./configure --prefix=$HOME/Downloads/libxslt-1.1.29_install
  make -j && make install
  ln -s $HOME/Downloads/libxslt-1.1.29_install/bin/xsltproc $HOME/bin
  popd
fi
#asciidoc
wget https://github.com/asciidoc/asciidoc/archive/8.6.9.zip
unzip 8.6.9.zip 
pushd asciidoc-8.6.9/
autoconf
./configure --prefix=$HOME/Downloads/asciidoc-8.6.9_install
make && make install
ln -sf $PWD/asciidoc.py ~/bin/asciidoc
ln -sf $PWD/a2x.py ~/bin/a2x
popd
#xmlto
wget http://pkgs.fedoraproject.org/repo/pkgs/xmlto/xmlto-0.0.28.tar.bz2/93bab48d446c826399d130d959fe676f/xmlto-0.0.28.tar.bz2
tar jxvf xmlto-0.0.28.tar.bz2 
pushd xmlto-0.0.28/
./configure --prefix=$HOME/Downloads/xmlto-0.0.28_install
make && make install
chmod +x $HOME/Downloads/xmlto-0.0.28_install/bin/xmlto
ln -sf $HOME/Downloads/xmlto-0.0.28_install/bin/xmlto ~/bin/
popd
#git
wget https://github.com/git/git/archive/v2.11.0.zip
unzip v2.11.0.zip
pushd git-2.11.0/
make configure && \
./configure --prefix=$HOME/Downloads/git-2.11.0_install --with-curl && \
make all doc && \
make install
ln -sf $HOME/Downloads/git-2.11.0_install/bin/git ~/bin
popd
#git bash autocompletion
echo 'source $HOME/Downloads/git-2.11.0/contrib/completion/git-completion.bash' >> $HOME/.bashrc


#####R#####
#libcurl
wget https://curl.haxx.se/download/curl-7.53.1.tar.gz
tar zxvf curl-7.53.1.tar.gz 
pushd curl-7.53.1/
./configure --prefix=$HOME/Downloads/curl-7.53.1_install && make -j && make install

#xz, lzma libs
wget http://tukaani.org/xz/xz-5.2.3.tar.gz
tar zxvf xz-5.2.3.tar.gz 
pushd xz-5.2.3/
./configure --prefix=$HOME/Downloads/xz-5.2.3_install && make -j && make install
popd
#R
wget https://cloud.r-project.org/src/base/R-3/R-3.3.2.tar.gz
tar zxvf R-3.3.2.tar.gz 
pushd R-3.3.2/
export CFLAGS="-I$HOME/Downloads/xz-5.2.3_install/include -I$HOME/Downloads/curl-7.53.1_install/include" &&\
export LDFLAGS="-L$HOME/Downloads/xz-5.2.3_install/lib -L$HOME/Downloads/curl-7.53.1_install/lib " &&\
./configure --disable-openmp --prefix=$HOME/Downloads/R-3.3.2_install &&\
make -j && make install
ln -s $HOME/Downloads/R-3.3.2_install/bin/* $HOME/bin
popd


#octave
curl https://ftp.gnu.org/gnu/octave/octave-4.2.1.tar.gz | tar zxv
pushd octave-4.2.1/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/isilon/Analysis/datainsights/users/guoy28/Downloads/OpenBLAS-0.2.19_install/lib && ./configure --prefix=$HOME/Downloads/octave-4.2.1_install --with-blas=$HOME/Downloads/OpenBLAS-0.2.19_install/lib --disable-readline
make -j && make install
ln -s $HOME/Downloads/octave-4.2.1_install/bin/octave $HOME/bin/
popd

#datamash
wget -O - http://ftp.gnu.org/gnu/datamash/datamash-1.1.0.tar.gz | tar xzf -
pushd datamash-1.1.0
./configure --prefix=$HOME/Downloads/datamash-1.1.0_install
make -j && make check && make install
ln -sf $HOME/Downloads/datamash-1.1.0_install/bin/datamash $HOME/bin
popd
echo 'export LIBRARY_PATH=$LD_LIBRARY_PATH:$LIBRARY_PATH' >> $HOME/.bashrc


#GNU parallel
curl https://ftp.gnu.org/gnu/parallel/parallel-20170522.tar.bz2 | tar jxvf -
pushd parallel-20170522
./configure --prefix=$HOME/Downloads/parallel-20170522_install
ln -sf $HOME/Downloads/parallel-20170522_install/bin/parallel ~/bin
popd

#bzip
wget -O- http://www.bzip.org/1.0.6/bzip2-1.0.6.tar.gz | tar zxvf -
pushd bzip2-1.0.6/
make
make install PREFIX=$HOME/Downloads/bzip2-1.0.6_install
popd
