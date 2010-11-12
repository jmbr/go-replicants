cmake_minimum_required(VERSION 2.6)

project(molecular-simulator C)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeScripts)

find_package(GSL REQUIRED)
find_package(Gnuplot REQUIRED)

add_definitions(-W -Wall -Wconversion -Wmissing-prototypes)
add_definitions(-Wstrict-prototypes -Wshadow -pedantic -std=c99)
add_definitions(-ggdb -O2)
add_definitions(-DHAVE_INLINE -D_XOPEN_SOURCE -D_BSD_SOURCE)
add_definitions(-DGNUPLOT_EXECUTABLE="${GNUPLOT_EXECUTABLE}")
include_directories(${GSL_INCLUDE_DIRS})

# This is an ugly hack to replace the CBLAS.  At some point,
# FindGSL.cmake ought to be rewritten.
list(REMOVE_ITEM GSL_LIBRARIES -lgslcblas)
list(APPEND GSL_LIBRARIES -lblas) # Use default BLAS (ATLAS).
link_directories(${GSL_LINK_DIRECTORIES})

add_executable(test-metropolis test-metropolis.c metropolis.c metropolis.h)
add_executable(test-protein test-protein.c protein.c protein.h utils.c utils.h)
add_executable(test-contact-map test-contact-map.c contact-map.c contact-map.h protein.c protein.h utils.c utils.h potential.c potential.h)

target_link_libraries(test-metropolis ${GSL_LIBRARIES})
target_link_libraries(test-protein ${GSL_LIBRARIES})
target_link_libraries(test-contact-map ${GSL_LIBRARIES})

enable_testing()
add_test(metropolis test-metropolis)
add_test(protein test-protein)
add_test(contact-map test-contact-map)
