# ------------------------------------------------------------------------------
# Project: Cadence Math Library
# Copyright (c) 2025, Onur Tuncer, PhD, Istanbul Technical University
#
# SPDX-License-Identifier: MIT
# License-Filename: LICENSE
# ------------------------------------------------------------------------------

find_package(Catch2 REQUIRED)

add_executable(${PROJECT_NAME}Tests
  test_Matrix.cpp
  test_Vector.cpp 
  test_Transform.cpp
)

target_include_directories(${PROJECT_NAME}Tests
  PRIVATE
    ${PROJECT_SOURCE_DIR}/  
)

target_link_libraries(${PROJECT_NAME}Tests PRIVATE
  Catch2::Catch2WithMain
)

enable_testing()
add_test(NAME ${PROJECT_NAME}Tests COMMAND ${PROJECT_NAME}Tests)
