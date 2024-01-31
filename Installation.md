This project uses a standard cmake build pattern to compile and run. Once the requirements (below) are set on your system, one should be able to, e.g., cd build; cmake ..; make

# Requirements

## CMAKE 3.14 minimum

## BOOST 1.66 minimum

tested on 1.84.0:
wget https://boostorg.jfrog.io/artifactory/main/release/1.84.0/source/boost_1_84_0.tar.gz
tar xf boost_1_84_0.tar.gz
cd boost_1_84_0
sudo ./boostrap.sh
sudo ./b2 install

## CGAL at least version 5.6 

wget https://github.com/CGAL/cgal/releases/download/v5.6/CGAL-5.6.tar.xz
tar xf CGAL-5.6.tar.xz
cd CGAL-5.6
cmake .
make install

## MPI

The code was tested with with openMPI 3.1 (anything later should also work)

## netcdf-cxx

Some of the database classes use netcdf to save data in a compressed form. Follow 
standard instructions to [install netcdf](https://docs.unidata.ucar.edu/nug/current/getting_and_building_netcdf.html), which will itself require zlib and hdf5. After this, one needs to [install netcdf-cxx](https://github.com/Unidata/netcdf-cxx4)

# Other people's code

Note that the include directory has files from the templatized c++ command line parser library (tclap), available from https://tclap.sourceforge.net/
