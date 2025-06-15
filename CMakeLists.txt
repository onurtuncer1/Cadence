
cmake_minimum_required(VERSION 3.26)
project(Cadence LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Header-only library
add_library(Cadence INTERFACE)
target_include_directories(Cadence INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)