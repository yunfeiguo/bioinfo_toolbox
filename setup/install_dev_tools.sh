mkdir -p $HOME/Downloads
cd $HOME/Downloads
NUM_THREADS=$( cat /proc/cpuinfo | grep -c '^processor')
USABLE_THREADS=$(expr $NUM_THREADS - 1)
#gcc, may a a few hours
wget http://www.netgull.com/gcc/releases/gcc-6.3.0/gcc-6.3.0.tar.gz
tar zxf gcc-6.3.0.tar.gz
pushd gcc-6.3.0/
./contrib/download_prerequisites
popd
mkdir gcc-6.3.0-objdir
pushd gcc-6.3.0-objdir/
../gcc-6.3.0/configure --prefix=$HOME/Downloads/gcc-6.3.0_install --enable-languages=c,c++,fortran,go --disable-multilib
make -j $USABLE_THREADS && make install
popd
echo 'export LD_LIBRARY_PATH=$HOME/Downloads/gcc-6.3.0_install/lib:$LD_LIBRARY_PATH' >> $HOME/.bashrc
echo 'export LD_LIBRARY_PATH=$HOME/Downloads/gcc-6.3.0_install/lib64:$LD_LIBRARY_PATH' >> $HOME/.bashrc
#cmake
#zlib
#python +setuptools+pip
wget https://www.python.org/ftp/python/2.7.12/Python-2.7.12.tgz
tar zxvf Python-2.7.12.tgz 
pushd Python-2.7.12/
./configure  --prefix=$HOME/Downloads/python-2.7.12_install
make -j 6 && make install
ln -sf $HOME/Downloads/python-2.7.12_install/bin/python $HOME/bin/

wget https://github.com/pypa/setuptools/archive/v29.0.1.zip
unzip v29.0.1.zip 
pushd setuptools-29.0.1/
python bootstrap.py 
python setup.py install
popd

wget https://pypi.python.org/packages/11/b6/abcb525026a4be042b486df43905d6893fb04f05aac21c32c638e939e447/pip-9.0.1.tar.gz#md5=35f01da33009719497f01a4ba69d63c9
tar zxvf pip-9.0.1.tar.gz 
pushd pip-9.0.1/
python setup.py install
ln -sf $HOME/Downloads/python-2.7.12_install/bin/pip ~/bin/
popd

popd
#perl
#boost
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
wget https://fedorahosted.org/releases/x/m/xmlto/xmlto-0.0.28.tar.bz2
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
make configure
./configure --prefix=$HOME/Downloads/git-2.11.0_install
make all doc
make install install-doc install-html
ln -sf $HOME/Downloads/git-2.11.0_install/bin/git ~/bin
popd
#git bash autocompletion
wget https://github.com/git/git/blob/master/contrib/completion/git-completion.bash
echo 'source $HOME/Downloads/git-completion.bash' >> $HOME/.bashrc

