if(${CMAKE_BUILD_TYPE} MATCHES "Debug")
    add_definitions(-DDEBUGFLAGUP)
endif(${CMAKE_BUILD_TYPE} MATCHES "Debug")

include_directories(incDirs)

add_library(simulation
    basicSimulation.cpp
    simulation.cpp
    mpiSimulation.cpp
    )
