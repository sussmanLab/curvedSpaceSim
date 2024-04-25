This project uses a standard cmake build pattern to compile and run. Once the requirements (below) are set on your system, one should be able to, e.g., cd build; cmake ..; make

# Sample compilation from a clean install of Ubuntu 22.04

Starting from a fresh copy of Ubuntu 22.04 as an example, the main branch  can be compiled by first installing Cmake, boost, and CGAL. The following commands will  do this:

##  CMake, Boost and CGAL:

sudo  apt install  cmake
wget https://boostorg.jfrog.io/artifactory/main/release/1.84.0/source/boost_1_84_0.tar.gz
tar xf boost_1_84_0.tar.gz
cd boost_1_84_0
sudo ./bootstrap.sh
sudo ./b2 install
cd ..
wget https://github.com/CGAL/cgal/releases/download/v5.6/CGAL-5.6.tar.xz
tar xf CGAL-5.6.tar.xz
cd CGAL-5.6
cmake .
make install
cd ..

## Other required packages:

Other dependencies can be install via  apt-get. We need netcdf-cxx (the legacy version,  4.2.X), which itself requires things like zlib, hdf5, and netcdf),  and the CGAL headers will need  gmp and mpfr. An MPI package  is needed as well. The following commands will do the trick:

sudo add-apt-repository universe
sudo apt-get update
sudo apt-get install zlib1g-dev libhdf5-dev libnetcdf-dev  netcdf-bin libnetcdf-cxx-legacy-dev libgmp-dev libmpfr-dev mpich

# Other people's code

Note that the include directory has files from the templatized c++ command line parser library (tclap), available from https://tclap.sourceforge.net/
