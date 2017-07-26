ROOT_DIR=$HOME/Downloads
pushd $ROOT_DIR
wget -O- http://www.graphviz.org/pub/graphviz/stable/SOURCES/graphviz-2.40.1.tar.gz | tar zxvf -
pushd graphviz-2.40.1
./configure --prefix=$ROOT_DIR/graphviz-2.40.1_install
make -j && make install
popd
ln -s $ROOT_DIR/graphviz-2.40.1_install/bin/* $HOME/bin/
