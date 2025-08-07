# ------------------------------------------------------------------------------
# Project: Cadence Math Library
# Copyright (c) 2025, Onur Tuncer, PhD, Istanbul Technical University
#
# SPDX-License-Identifier: MIT
# License-Filename: LICENSE
# ------------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.22)

project(Cadence
  VERSION 1.0.0
  LANGUAGES CXX
  DESCRIPTION "Header-only math library for C++"
)

# Use cmake/ folder for custom modules (if any)
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Create header-only interface library
add_library(${PROJECT_NAME} INTERFACE)
add_library(Cadence::Cadence ALIAS ${PROJECT_NAME})

# Public include path for the headers
target_include_directories(${PROJECT_NAME}
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
)

# Optional: Build docs, run analysis, add tests
if (BUILD_DOCS)
  add_subdirectory(doc)
endif()

if (ENABLE_STATIC_ANALYSIS)
  include(cppcheck)
endif()

if (BUILD_TESTS)
  enable_testing()
  add_subdirectory(test)
endif()
