if(${CMAKE_BUILD_TYPE} MATCHES "Debug")
    add_definitions(-DDEBUGFLAGUP)
endif(${CMAKE_BUILD_TYPE} MATCHES "Debug")

include_directories(incDirs)

add_library(updaters
    baseUpdater.cpp
    gradientDescent.cpp
    fireMinimization.cpp
    velocityVerletNVE.cpp
    noseHooverNVT.cpp
    diffusiveSPP.cpp
    )
