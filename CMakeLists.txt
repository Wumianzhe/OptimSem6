cmake_minimum_required(VERSION 3.12)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(metopt)
set(CMAKE_BUILD_TYPE Debug)
add_subdirectory(libmatrix)
add_executable(lab1 main.cpp)
target_link_libraries(lab1 simplemat)
