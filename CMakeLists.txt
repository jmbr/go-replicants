cmake_minimum_required(VERSION 2.6)

project(molecular-simulator C)

enable_language(C CXX)          # C++ is needed for finding OpenMP.

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeScripts)

find_package(GSL REQUIRED)
find_package(Gnuplot REQUIRED)
find_package(OpenMP REQUIRED)

set(SOURCES
    protein.c protein.h
    contact-map.c contact-map.h
    utils.c utils.h
    potential.c potential.h
    geometry.c geometry.h
    simulation.c simulation.h)

add_definitions(-W -Wall -Wconversion -Wmissing-prototypes)
add_definitions(-Wstrict-prototypes -Wshadow -pedantic -std=c99)
add_definitions(-O2 -ggdb -DHAVE_INLINE -DNDEBUG -DDEBUG_LEVEL=0)
# add_definitions(-ggdb -O0 -DDEBUG_LEVEL=1)
add_definitions(-D_XOPEN_SOURCE -D_BSD_SOURCE)
add_definitions(-DGNUPLOT_EXECUTABLE="${GNUPLOT_EXECUTABLE}")

include_directories(${GSL_INCLUDE_DIRS})

# This is an ugly hack to replace the CBLAS.  At some point,
# FindGSL.cmake ought to be rewritten.
list(REMOVE_ITEM GSL_LIBRARIES -lgslcblas)
list(APPEND GSL_LIBRARIES -lblas) # Use default BLAS (ATLAS).
list(APPEND GSL_LIBRARIES -pg)
link_directories(${GSL_LINK_DIRECTORIES})

add_executable(test-protein test-protein.c ${SOURCES})
add_executable(test-contact-map test-contact-map.c ${SOURCES})
add_executable(test-simulation test-simulation.c ${SOURCES})

target_link_libraries(test-protein ${GSL_LIBRARIES})
target_link_libraries(test-contact-map ${GSL_LIBRARIES})
target_link_libraries(test-simulation ${GSL_LIBRARIES})

enable_testing()
add_test(protein test-protein)
add_test(contact-map test-contact-map)
add_test(simulation test-simulation)
