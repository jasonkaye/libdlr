# set minimum cmake version
cmake_minimum_required(VERSION 3.12)

# project name and language
project(id_dist LANGUAGES Fortran)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

file(GLOB id_SRC CONFIGURE_DEPENDS "./src/*.f")

add_library(id_dist ${id_SRC})
target_link_libraries(id_dist PRIVATE BLAS::BLAS LAPACK::LAPACK)
install(TARGETS id_dist LIBRARY)
install(TARGETS id_dist ARCHIVE)
