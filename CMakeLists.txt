cmake_minimum_required(VERSION 3.16)

project(GaussDistributionTraceGeneration LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(gauss_trace
    src/curves.cpp
    src/simulation.cpp
)

target_include_directories(gauss_trace
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

add_executable(generate generate.cpp)
target_link_libraries(generate PRIVATE gauss_trace)

add_executable(smoke_tests tests/smoke_tests.cpp)
target_link_libraries(smoke_tests PRIVATE gauss_trace)

enable_testing()
add_test(NAME smoke_tests COMMAND smoke_tests)
