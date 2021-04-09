# set minimum cmake version
cmake_minimum_required(VERSION 3.12)

# project name and language
project(id_dist LANGUAGES Fortran)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

add_subdirectory(./id_dist)

file(GLOB dlr_SRC CONFIGURE_DEPENDS "./src/*.f*" "./src/utils_external/*.f*")

add_library(dlr SHARED ${dlr_SRC})
target_link_libraries(dlr PUBLIC id_dist BLAS::BLAS LAPACK::LAPACK)

add_executable(dlr_mf_sc_test ./test/dlr_mf_sc_test.f90)
target_link_libraries(dlr_mf_sc_test PRIVATE dlr)

add_executable(dlr_it_sc_test ./test/dlr_it_sc_test.f90)
target_link_libraries(dlr_it_sc_test PRIVATE dlr)
