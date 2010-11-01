cmake_minimum_required(VERSION 2.6)

project(molecular-simulator)

include("FindGSL.cmake")

enable_language(C)

add_definitions(-Wall -pedantic -std=c99)
include_directories(${INCLUDE_DIRECTORIES} ${GSL_INCLUDE_DIRS})
link_directories(${LINK_DIRECTORIES} ${GSL_LINK_DIRECTORIES})
add_executable(test-metropolis test-metropolis.c metropolis.c metropolis.h)
target_link_libraries(test-metropolis ${GSL_LIBRARIES})

enable_testing()
add_test(metropolis test-metropolis)
