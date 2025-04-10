cmake_minimum_required(VERSION 3.11.0)

project(curvedSpaceSimulation LANGUAGES CXX)

#note: CGAL (which I will get rid of someday) needs at least c++14; some advanced gcc compilers will throw errors if you don't set the cpp standard to 17, though
set(CMAKE_CXX_STANDARD 17)

#Release by default
set(CMAKE_BUILD_TYPE "Release")
#set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
#set(CMAKE_VERBOSE_MAKEFILE on)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

# If on MacOS, assume that pthreads in built-in
IF(APPLE)
    set(CMAKE_THREAD_LIBS_INIT "-lpthread")
    set(CMAKE_HAVE_THREADS_LIBRARY 1)
    set(CMAKE_USE_WIN32_THREADS_INIT 0)
    set(CMAKE_USE_PTHREADS_INIT 1)
    set(THREADS_PREFER_PTHREAD_FLAG ON)
ENDIF()


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_MATH_DISABLE_FLOAT128 -fopenmp -frounding-math -Wall -Wshadow -Wno-unused-function -Wno-sign-compare")

if(${CMAKE_BUILD_TYPE} MATCHES "Debug")
    add_definitions(-DDEBUGFLAGUP)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g")
else()
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -pg") #for gprof profiling
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ffast-math -march=native")
    add_definitions(-DNDEBUG)
endif()
message("cxx flags = " "${CMAKE_CXX_FLAGS}")

find_package(CGAL REQUIRED)
find_package(Boost REQUIRED)
find_package(HDF5 REQUIRED)
find_package(MPI REQUIRED)

include_directories(
            ${CMAKE_SOURCE_DIR}/inc
            ${CMAKE_SOURCE_DIR}/src
            ${CMAKE_SOURCE_DIR}/src/utility
            ${CMAKE_SOURCE_DIR}/src/forces
            ${CMAKE_SOURCE_DIR}/src/updaters
            ${CMAKE_SOURCE_DIR}/src/simulation
            ${CMAKE_SOURCE_DIR}/src/models
            ${CMAKE_SOURCE_DIR}/src/analysis
            ${CMAKE_SOURCE_DIR}/src/databases
            ${MPI_INCLUDE_PATH}
   )

add_subdirectory(src)

set(myLibs
    model
    updaters
    databases
    simulation
    utility
    forces
    )

link_directories($ENV{LD_LIBRARY_PATH})

set(SOURCES
    )
set(HEADERS
    )

# list the names of cpp files corresponding to linked executables you'd like...
foreach(ARG
        meshTesting
        flatSpaceSimulation
        curvedSpaceSimulation
        multirankSimulation
        curvedSpaceNVTSim
        tangentialSimulation
        absorbingSimulation
	)

add_executable("${ARG}.out" "${ARG}.cpp" )
target_link_libraries("${ARG}.out"
    ${myLibs}
    #OpenMP::OpenMP_CXX
    CGAL::CGAL
    ${MPI_LIBRARIES}
    )
endforeach()
