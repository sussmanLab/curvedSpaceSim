if(${CMAKE_BUILD_TYPE} MATCHES "Debug")
    add_definitions(-DDEBUGFLAGUP)
endif(${CMAKE_BUILD_TYPE} MATCHES "Debug")

include_directories(incDirs)

add_library(databases
    DatabaseNetCDF.cpp
    vectorValueDatabase.cpp
    simpleModelDatabase.cpp
    )
