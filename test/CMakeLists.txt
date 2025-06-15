

# Fetch Catch2
include(FetchContent)
FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v3.5.3 # or latest stable
)
FetchContent_MakeAvailable(Catch2)

# Tests
enable_testing()
add_executable(vec_tests test/test_vector.cpp)
target_link_libraries(vector_tests PRIVATE Cadence Catch2::Catch2WithMain)

include(CTest)
catch_discover_tests(vector_tests)