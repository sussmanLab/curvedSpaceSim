if(${CMAKE_BUILD_TYPE} MATCHES "Debug")
    add_definitions(-DDEBUGFLAGUP)
endif(${CMAKE_BUILD_TYPE} MATCHES "Debug")

include_directories(incDirs)

add_library(databases
    baseHDF5Database.cpp
    vectorValueDatabase.cpp
    simpleModelDatabase.cpp
    )

target_include_directories(databases PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_libraries(databases PUBLIC ${HDF5_LIBRARIES})
