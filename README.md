# curvedSpaceSimulation

Lorum ipsum, to the max

# Requirements

# CMAKE 3.11 minimum

# BOOST 1.66 minimum

tested on 1.84.0:
wget https://boostorg.jfrog.io/artifactory/main/release/1.84.0/source/boost_1_84_0.tar.gz
tar xf boost_1_84_0.tar.gz
cd boost_1_84_0
sudo ./boostrap.sh
sudo ./b2 install

# CGAL at least version 5.6 

wget https://github.com/CGAL/cgal/releases/download/v5.6/CGAL-5.6.tar.xz
tar xf CGAL-5.6.tar.xz
cd CGAL-5.6
cmake .
make install

# MPI
Tested with openMPI 3.1 (anything ater should also work)
