cmake_minimum_required(VERSION 3.3.0)
project(piap)

set(EXECUTABLE_NAME "piap")

include_directories(src)
file(GLOB_RECURSE SOURCES "src/*.cpp")

add_executable(${EXECUTABLE_NAME} ${SOURCES})

set_property(TARGET ${EXECUTABLE_NAME} PROPERTY CXX_STANDARD 14)
set_property(TARGET ${EXECUTABLE_NAME} PROPERTY CXX_STANDARD_REQUIRED ON)
