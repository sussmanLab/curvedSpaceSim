if(${CMAKE_BUILD_TYPE} MATCHES "Debug")
    add_definitions(-DDEBUGFLAGUP)
endif(${CMAKE_BUILD_TYPE} MATCHES "Debug")

include_directories(incDirs)

add_library(utility
    functionUtilities.cpp
    meshUtilities.cpp
    hyperRectangularCellList.cpp
    baseNeighborStructure.cpp
    cellListNeighborStructure.cpp
    noiseSource.cpp
    )
