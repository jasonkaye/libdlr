# set minimum cmake version
cmake_minimum_required(VERSION 3.12)

# project name and language
project(dlr LANGUAGES Fortran)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ../lib)
add_subdirectory(./id_dist)

# Fortran library: dlr

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY lib)

file(GLOB dlr_SRC CONFIGURE_DEPENDS "./src/*.f*" "./src/utils_external/*.f*")

add_library(dlr SHARED ${dlr_SRC})
target_link_libraries(dlr PUBLIC id_dist BLAS::BLAS LAPACK::LAPACK)
install(TARGETS dlr LIBRARY)

# Test programs

enable_testing()

add_executable(dlr_mf_sc_test ./test/dlr_mf_sc_test.f90)
target_link_libraries(dlr_mf_sc_test PRIVATE dlr)
add_test(NAME dlr_mf_sc_test COMMAND dlr_mf_sc_test)

add_executable(dlr_it_sc_test ./test/dlr_it_sc_test.f90)
target_link_libraries(dlr_it_sc_test PRIVATE dlr)
add_test(NAME dlr_it_sc_test COMMAND dlr_it_sc_test)

# C-library: dlr_c

file(GLOB dlr_c_SRC CONFIGURE_DEPENDS "./src/dlr_c/*.f*")
add_library(dlr_c SHARED ${dlr_c_SRC})
target_link_libraries(dlr_c PUBLIC dlr id_dist BLAS::BLAS LAPACK::LAPACK)
install(TARGETS dlr_c LIBRARY)
set_target_properties(dlr_c PROPERTIES INSTALL_RPATH "@loader_path/")

# Python module: pydlr

find_package(Python COMPONENTS Interpreter Development)
FILE(GLOB_RECURSE pydlr_SRC RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.py)
foreach(f ${pydlr_SRC})
 configure_file(${f} ${f} COPYONLY)
endforeach()
install(FILES pydlr/pydlr.py DESTINATION pydlr)

add_test(NAME python-test_kernel
  COMMAND ${Python_EXECUTABLE} -m nose ./test/test_kernel.py -v
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pydlr
  )
