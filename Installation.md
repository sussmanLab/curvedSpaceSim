This project uses a standard cmake build pattern to compile and run. Once the requirements (below) are set on your system, one should be able to, e.g., cd build; cmake ..; make

# Sample compilation from a clean install of Ubuntu 22.04

Starting from a fresh copy of Ubuntu 22.04 as an example, the main branch  can be compiled by first installing Cmake, boost, and CGAL. Start in a directory that you don't mind adding some tar files and software folders to, and then run the following commands:

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

## Compiling the curvedSpaceSimulation executables

With all of the dependencies installed, compilation of the code is done with a standard out-of-source build: From the min repository directory: 

cd build

cmake ..

make

Note that depending on how your local system is set up, you may need to give cmake additional hints about where various libraries are.
If you want to write a new executable cpp file, of course just add it to the main CMakeLists.txt file in the "foreach(ARG " list around line 70

# Installing dependencies on a generic linux system (without admin privileges)

If you are unable to sudo, the following will walk through the compilation and installation of the required packages from the command line.

## setting up a local library, include, and bin location, and pointing relevant path variables there:

For convenience, you might consider installing all of these dependencies in a shared local directory, and then updating your various path variables to point to it (making CMake compilation a bit easier)...replace zshrc with bashrc as needed

mkdir -p $HOME/.local

mkdir -p $HOME/.local/bin

mkdir -p $HOME/.local/include

mkdir -p $HOME/.local/lib

echo 'export PATH="$PATH:$HOME/.local/bin"' >> ~/.zshrc

echo 'export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/.local/lib"'  >> ~/.zshrc

echo 'export LIBRARY_PATH="$LIBRARY_PATH:$HOME/.local/lib"' >> ~/.zshrc

echo 'export CPATH="$CPATH:$HOME/.local/include"' >> ~/.zshrc

echo 'export GMP_LIBRARIES=$HOME/.local/lib'  >> ~/.zshrc

echo 'export GMP_INCLUDE_DIR=$HOME/.local/include' >> ~/.zshrc

echo 'export MPFR_LIBRARIES=$HOME/.local/lib' >> ~/.zshrc

echo 'export MPFR_INCLUDE_DIR=$HOME/.local/include'  >> ~/.zshrc

then do a ``source ~/.zshrc''. CMake may still have some problems finding all of these; I'm not enough of a CMake wizard to know what the problem might be (and it might depend on the version of cmake you're using). The upshot is that you might have to monkey a bit with the CMakeCache.txt file to point cmake to the right folders and libraries. 

Anyway, we now install all of the dependencies. Go to some root directory where you don't mind accumulating tar files and software directories.


## install boost

curl -L -O https://boostorg.jfrog.io/artifactory/main/release/1.84.0/source/boost_1_84_0.tar.gz

tar xf boost_1_84_0.tar.gz

cd boost_1_84_0

./bootstrap.sh --prefix=$HOME/.local

./b2 install

cd ..


## install gmp

curl -L -O https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz

tar xf gmp-6.3.0.tar.xy

cd gmp-6.3.0

./configure --prefix=$HOME/.local

make

make check

make install

cd ..


## install MPFR

curl -L -O  "https://www.mpfr.org/mpfr-current/mpfr-4.2.1.tar.xz"

tar -xvf mpfr-4.2.1.tar.xz

cd mpfr-4.2.1

./configure --prefix=$HOME/.local

make

make check

make install

cd ..


## install CGAL

curl -L -O https://github.com/CGAL/cgal/releases/download/v5.6.1/CGAL-5.6.1.tar.xz

tar -xvf CGAL-5.6.1.tar.xz 

cd CGAL-5.6.1

cmake -DCMAKE_INSTALL_PREFIX=$HOME/.local -DCMAKE_BUILD_TYPE=Release .

make

make install

## install open-MPI

curl -L -O https://download.open-mpi.org/release/open-mpi/v5.0/openmpi-5.0.3.tar.gz

tar -xvf openmpi-5.0.3.tar.gz   

cd openmpi-5.0.3

./configure --prefix=$HOME/.local

make all

make install

## install zlib

curl -L -O  "http://www.zlib.net/zlib-1.3.1.tar.gz"

tar axf zlib-1.3.1.tar.gz

cd zlib-1.3.1

./configure --prefix=$HOME/.local

make

make check

make install

cd ..

## install hdf5

curl -L -O  "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.7/src/hdf5-1.10.7.tar.gz"

tar xf hdf5-1.10.5.tar.gz

cd hdf5-1.10.5

./configure --prefix=$HOME/.local --enable-cxx

make

make check

make install

## install netcdf

curl -L -O  "https://downloads.unidata.ucar.edu/netcdf-c/4.9.2/netcdf-c-4.9.2.tar.gz"

tar axf "netcdf-c-4.9.2.tar.gz"

cd netcdf-c-4.9.2

./configure --prefix=$HOME/.local --enable-netcdf-4

make

make check

make install

cd ..

## install netcdf-cxx

curl -L -O "https://downloads.unidata.ucar.edu/netcdf-cxx/4.2/netcdf-cxx-4.2.tar.gz"

tar -xf netcdf-cxx-4.2.tar.gz

cd netcdf-cxx-4.2

./configure --prefix=$HOME/.local

make

make check

make install

cd ..

# Installing on MacOS

The following was tested on a fresh account of macOS 14.3.1 running on an M1 chip. XCode had already been set up, with the  Apple clang version 15.0.0 toolchain. If you have access to homebrew many of the dependencies can be "brew installed" -- google will be your friend here (unfortunately, the current developers do not have ready access to an apple device on which they can check whether the brew methods all work together. So, if you,too,  want to compile dependencies from the command line, read on.

## local CMake install

On a mac the easiest thing is just to "brew install cmake"... if you do not have access to homebrew, or want to do everything without administrator abilities on your mac, you can do the following (with your user name in place of "userName") instead:

$ curl --silent --location --retry 3 "https://github.com/Kitware/CMake/releases/download/v3.28.3/cmake-3.28.3-macos-universal.dmg" --output ~/Downloads/CMake/cmake-macos.dmg
$ yes | PAGER=cat hdiutil attach -quiet -mountpoint /Volumes/cmake-macos ~/Downloads/CMake/cmake-macos.dmg
$ mkdir /Users/userName/Applications
$ cp -R /Volumes/cmake-macos/CMake.app /Users/userName/Applications/
$ hdiutil detach /Volumes/cmake-macos

Then just make a "cmake" alias in your .bashrc or .zshrc file that points to the bin directory of the CMake/Contents/ directory in that local applications folder

## Remaining steps

For the most part you should be able to follow the same steps as the generic linux install in the previous section, with two issues

### netcdf-cxx

Compiling on an Apple device with clang requires some modifications to the steps above. Start out as before by downloading and going into the directory:

curl -L -O "https://downloads.unidata.ucar.edu/netcdf-cxx/4.2/netcdf-cxx-4.2.tar.gz"

tar -xf netcdf-cxx-4.2.tar.gz

cd netcdf-cxx-4.2

Now, before configuring, use a text editor to go into the "configure" file and find the line that says "# Create the VERSION file, which contains the package version from" (this should be line 2330). Below that is an uncommented line that says "echo -n 4.2>VERSION"; change this to "echo -n 4.2>VERSION.txt". It's janky, and I don't know why that matters, but it works. Once that change is made, procede as before:

./configure --prefix=$HOME/.local

make

make check

make install

cd ..

### open-MPI

The instructions for open-MPI installation above will work just fine. There are, however, some known issues with using clang to compile open-mpi (leading to problems in which cmake will complain about not being able to find MPI_CXX). You may be able to use some CMake magic to get around this. Alternatives include (a) downloading and compiling open-mpi with gcc rather than clang (?), or (b) not using MPI on your mac. 

To pursue option (b), a few changes in the root CMakeLists.txt file need to be made:

(1) On line 23 get rid of the "-fopenmp" flag

(2) Comment out line 37 (the find_pacakge(MPI) line)

(3) Comment out (or delete) line 74 (so that multirankSimulation doesn't try to compile)

(4) Comment out (or delete) line 83 (so that ${MPI_LIBRARIES} aren't linked)

# Other people's code

Note that the include directory has files from the templatized c++ command line parser library (tclap), available from https://tclap.sourceforge.net/
