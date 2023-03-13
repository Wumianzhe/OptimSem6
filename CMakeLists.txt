cmake_minimum_required(VERSION 3.20)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# add_compile_options(-ffinite-math-only -fno-rounding-math)
add_compile_options($<$<CONFIG:Debug>:-fsanitize=address>)
add_link_options($<$<CONFIG:Debug>:-fsanitize=address>)

project(metopt)
set(CMAKE_BUILD_TYPE Debug)

include_directories("${PROJECT_SOURCE_DIR}/inc")

add_subdirectory(libmatrix)
add_executable(lab main.cpp)
target_link_libraries(lab simplemat)
