cmake_minimum_required(VERSION 2.6)

project(molecular-simulator C)

enable_language(C CXX Fortran)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMake)

find_package(GSL REQUIRED)
find_package(Gnuplot)
# find_package(BLAS)
find_package(OpenMP)
# find_package(MKL)
if(OPENMP_FOUND)
  list(APPEND CMAKE_C_FLAGS ${OpenMP_C_FLAGS})
endif()
# if(BLAS_FOUND)
#   list(REMOVE_ITEM GSL_LIBRARIES -lgslcblas)
# endif()

add_definitions(-D_XOPEN_SOURCE -D_BSD_SOURCE)
add_definitions(-W -Wall -Wconversion -Wmissing-prototypes)
add_definitions(-Wstrict-prototypes -Wshadow -pedantic -std=c99)
add_definitions(-ggdb -DHAVE_INLINE -DNDEBUG -DDEBUG_LEVEL=0)
add_definitions(-DGNUPLOT_EXECUTABLE="${GNUPLOT_EXECUTABLE}")

set(SIMULATOR_SOURCE_FILES protein.c protein.h movements.c movements.h
  contact-map.c contact-map.h utils.c utils.h potential.c
  potential.h geometry.c geometry.h simulation.c simulation.h
  replicas.c replicas.h)

add_library(simulator ${SIMULATOR_SOURCE_FILES})

add_executable(molecular-simulator molecular-simulator.c)
add_executable(molecular-viewer molecular-viewer.c)
add_executable(molecular-player molecular-player.c)
add_executable(eval-potential eval-potential.c)
add_executable(test-protein test-protein.c)
add_executable(test-contact-map test-contact-map.c)
add_executable(test-replicas test-replicas.c)

set(TARGETS molecular-simulator molecular-viewer molecular-player eval-potential)

set_target_properties(${TARGETS}
  PROPERTIES LINK_FLAGS "${GSL_LINKER_FLAGS} ${BLAS_LINKER_FLAGS}")

target_link_libraries(molecular-simulator simulator ${GSL_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(molecular-viewer simulator ${GSL_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(molecular-player simulator ${GSL_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(eval-potential simulator ${GSL_LIBRARIES} ${BLAS_LIBRARIES})

enable_testing()
add_test(protein test-protein)
add_test(contact-map test-contact-map)
add_test(replicas test-replicas)
set_tests_properties(protein contact-map replicas
  PROPERTIES LINK_FLAGS "${GSL_LINKER_FLAGS} ${BLAS_LINKER_FLAGS}")
target_link_libraries(test-protein simulator ${GSL_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(test-contact-map simulator ${GSL_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(test-replicas simulator ${GSL_LIBRARIES} ${BLAS_LIBRARIES})

include(CPack)
