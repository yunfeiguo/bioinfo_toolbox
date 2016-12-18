mkdir -p $HOME/Downloads
cd $HOME/Downloads
#gcc
#cmake
#zlib
#python
#perl
#boost
#asciidoc
wget https://github.com/asciidoc/asciidoc/archive/8.6.9.zip
unzip 8.6.9.zip 
pushd asciidoc-8.6.9/
autoconf
./configure --prefix=$HOME/Downloads/asciidoc-8.6.9_install
make && make install
ln -sf $HOME/Downloads/asciidoc-8.6.9_install/asciidoc.py ~/bin/asciidoc
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
