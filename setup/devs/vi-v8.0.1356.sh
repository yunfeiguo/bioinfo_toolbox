git clone https://github.com/vim/vim.git
pushd vim
git checkout v8.0.1356
pushd src/
./configure --prefix=$HOME/Downloads/vim_v8.0.1356_install --enable-pythoninterp=yes && make -j && make install && ln -sf $HOME/Downloads/vim_v8.0.1356_install/bin/* $HOME/bin
make uninstall && make clean && ./configure --prefix=$HOME/Downloads/vim_v8.0.1356_install --enable-pythoninterp=yes && make -j && make install && ln -sf $HOME/Downloads/vim_v8.0.1356_install/bin/* $HOME/bin
popd
popd
