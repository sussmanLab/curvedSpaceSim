if(${CMAKE_BUILD_TYPE} MATCHES "Debug")
    add_definitions(-DDEBUGFLAGUP)
endif(${CMAKE_BUILD_TYPE} MATCHES "Debug")

include_directories(incDirs)

add_library(model
    simpleModel.cpp
    mpiModel.cpp
    euclideanSpace.cpp
    triangulatedMeshSpace.cpp
    openMeshSpace.cpp
    absorbingOpenMeshSpace.cpp
    tangentialOpenMeshSpace.cpp
    )
