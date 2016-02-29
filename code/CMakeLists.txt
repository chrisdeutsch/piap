cmake_minimum_required(VERSION 3.3.0)
project(piap)

include_directories(src)

add_executable(coulomb2d src/coulomb2d.cpp)

set_property(TARGET coulomb2d PROPERTY CXX_STANDARD 14)
set_property(TARGET coulomb2d PROPERTY CXX_STANDARD_REQUIRED ON)
