pushd $HOME/Downloads
wget -O- http://www.cpan.org/src/5.0/perl-5.26.0.tar.gz | tar zxvf -
pushd perl-5.26.0
./Configure -des -Dprefix=$HOME/Downloads/perl-5.26.0_install && make -j && make install
popd
ln -sf $PWD/perl-5.26.0_install/bin/* $HOME/bin
popd
