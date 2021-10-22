# SPM_Project

On Linux or Mac OS X, another option is to symlink or copy the Eigen folder into /usr/local/include/. This way, you can compile the program with:

cp eigen-3.4.0 /usr/local/include
cd /usr/local/include
sudo ln -sf eigen-3.4.0/Eigen/ Eigen

or use -I eigen-3.4.0