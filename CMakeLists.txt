cmake_minimum_required(VERSION 2.6)

project(molecular-simulator C)

include(FindGSL.cmake)

add_definitions(-W -Wall -Wconversion -Wmissing-prototypes -Wstrict-prototypes -Wshadow -std=c99 -pedantic -ggdb -O2)
include_directories(${INCLUDE_DIRECTORIES} ${GSL_INCLUDE_DIRS})

link_directories(${LINK_DIRECTORIES} ${GSL_LINK_DIRECTORIES})

add_executable(test-metropolis test-metropolis.c metropolis.c metropolis.h)
add_executable(test-molecule test-molecule.c molecule.c molecule.h molecule-private.h xalloc.c xalloc.h error.c error.h)

target_link_libraries(test-metropolis ${GSL_LIBRARIES})

enable_testing()
add_test(metropolis test-metropolis)
add_test(molecule test-molecule)
