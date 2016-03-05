cmake_minimum_required(VERSION 3.3.0)
project(piap)

include_directories(src)

# Coulomb with hard core in 2 dimension
set(SOURCES_2D
    src/coulomb2d.cpp
    src/Particle.cpp)

add_executable(coulomb2d ${SOURCES_2D})

set_property(TARGET coulomb2d PROPERTY CXX_STANDARD 14)
set_property(TARGET coulomb2d PROPERTY CXX_STANDARD_REQUIRED ON)
