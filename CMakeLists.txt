cmake_minimum_required(VERSION 2.6)

project(molecular-simulator C)

enable_language(C CXX)          # C++ is needed for finding OpenMP.

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMake)

find_package(GSL REQUIRED)
find_package(Gnuplot)
find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif()

add_library(simulator protein.c protein.h movements.c movements.h
contact-map.c contact-map.h linspace.h utils.c utils.h potential.c
potential.h geometry.c geometry.h simulation.c simulation.h replicas.c
replicas.h)

add_definitions(-D_XOPEN_SOURCE -D_BSD_SOURCE)
add_definitions(-W -Wall -Wconversion -Wmissing-prototypes)
add_definitions(-Wstrict-prototypes -Wshadow -pedantic -std=c99)
add_definitions(-ggdb -O2 -DHAVE_INLINE -DNDEBUG -DDEBUG_LEVEL=0)
add_definitions(-DGNUPLOT_EXECUTABLE="${GNUPLOT_EXECUTABLE}")

include_directories(${GSL_INCLUDE_DIRS})

# This is an ugly hack to replace the CBLAS.  At some point,
# FindGSL.cmake ought to be rewritten.
list(REMOVE_ITEM GSL_LIBRARIES -lgslcblas)
list(APPEND GSL_LIBRARIES -lblas) # Use default BLAS (ATLAS).
list(APPEND GSL_LIBRARIES)
link_directories(${GSL_LINK_DIRECTORIES})

add_executable(molecular-simulator molecular-simulator.c)
add_executable(molecular-viewer molecular-viewer.c)
add_executable(molecular-player molecular-player.c)
add_executable(eval-potential eval-potential.c)
add_executable(test-protein test-protein.c)
add_executable(test-contact-map test-contact-map.c)
add_executable(test-replicas test-replicas.c)

# set_target_properties(simulator PROPERTIES LINK_FLAGS -pg)
# set_target_properties(molecular-simulator PROPERTIES LINK_FLAGS -pg)

target_link_libraries(molecular-simulator simulator ${GSL_LIBRARIES})
target_link_libraries(molecular-viewer simulator ${GSL_LIBRARIES})
target_link_libraries(molecular-player simulator ${GSL_LIBRARIES})
target_link_libraries(eval-potential simulator ${GSL_LIBRARIES})
target_link_libraries(test-protein simulator ${GSL_LIBRARIES})
target_link_libraries(test-contact-map simulator ${GSL_LIBRARIES})
target_link_libraries(test-replicas simulator ${GSL_LIBRARIES})

enable_testing()
add_test(protein test-protein)
add_test(contact-map test-contact-map)
add_test(replicas test-replicas)

include(CPack)
