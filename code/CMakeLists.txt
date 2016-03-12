cmake_minimum_required(VERSION 3.3.0)
project(piap)

include_directories(src)


# Coulomb with hard core in 2 dimensions
set(SOURCES_2D
    src/coulomb2d.cpp
    src/Particle.cpp
    src/Common.cpp)

add_executable(coulomb2d ${SOURCES_2D})


# Calculate observables
set(SOURCES_2D_OBS
    src/coulomb2d_obs.cpp
    src/Particle.cpp
    src/Common.cpp)

add_executable(coulomb2d_obs ${SOURCES_2D_OBS})


# Coulomb with hard core in 3 dimensions
set(SOURCES_3D
    src/coulomb3d.cpp
    src/Particle.cpp
    src/Common.cpp)

add_executable(coulomb3d ${SOURCES_3D})


set(TARGETS
    coulomb2d
    coulomb2d_obs
    coulomb3d)

set_property(TARGET ${TARGETS} PROPERTY CXX_STANDARD 14)
set_property(TARGET ${TARGETS} PROPERTY CXX_STANDARD_REQUIRED ON)
