# SPM_Project


mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX='../spectrainstall' -DCMAKE_PREFIX_PATH='../eigen-3.4.0' -DBUILD_TESTS=FALSE
make all && make tests && make install
