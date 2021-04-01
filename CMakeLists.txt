# set minimum cmake version
cmake_minimum_required(VERSION 3.12)

# project name and language
project(id_dist LANGUAGES Fortran)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

add_subdirectory(./id_dist)

file(GLOB dlr_SRC CONFIGURE_DEPENDS "./src/*.f*" "./src/utils_external/*.f*")

# Fortran library

add_library(dlr SHARED ${dlr_SRC})
target_link_libraries(dlr PUBLIC id_dist BLAS::BLAS LAPACK::LAPACK)

# C-library

file(GLOB dlr_c_SRC CONFIGURE_DEPENDS "./src/dlr_c/*.f*")
add_library(dlr_c SHARED ${dlr_c_SRC})
target_link_libraries(dlr_c PUBLIC dlr id_dist BLAS::BLAS LAPACK::LAPACK)

# Test programs

add_executable(dlr_mf_sc_test ./test/dlr_mf_sc_test.f90)
target_link_libraries(dlr_mf_sc_test PRIVATE dlr)

add_executable(dlr_it_sc_test ./test/dlr_it_sc_test.f90)
target_link_libraries(dlr_it_sc_test PRIVATE dlr)
