# set minimum cmake version
cmake_minimum_required(VERSION 3.12)

# project name and language
project(dlr LANGUAGES Fortran)

option(with_python "Install python module: pydlr" OFF)

include(CTest)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ../lib)
add_subdirectory(./id_dist)

# Fortran library: dlr

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY lib)

file(GLOB dlr_SRC CONFIGURE_DEPENDS "./src/*.f*" "./src/utils_external/*.f*")

add_library(dlr STATIC ${dlr_SRC})
target_link_libraries(dlr PUBLIC id_dist BLAS::BLAS LAPACK::LAPACK)
install(TARGETS dlr LIBRARY)

# Test programs

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

# Macro for finding Python modules

macro(find_python_module _module_name)
  execute_process(
    COMMAND
      ${Python3_EXECUTABLE} "-c" "import ${_module_name}; print(${_module_name}.__version__)"
    OUTPUT_VARIABLE _stdout
    ERROR_VARIABLE _stderr
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_STRIP_TRAILING_WHITESPACE
    )

  if(_stderr MATCHES "ModuleNotFoundError")
    message(SEND_ERROR "Python module \"${_module_name}\" not found")
  else()
    message(STATUS "Found Python module ${_module_name} v${_stdout}")
  endif()
endmacro()

# Python module: pydlr

if(with_python)

  find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
  find_python_module("numpy")
  find_python_module("cffi")
  find_python_module("nose")

  FILE(GLOB_RECURSE pydlr_SRC RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.py)
  foreach(f ${pydlr_SRC})
   configure_file(${f} ${f} COPYONLY)
  endforeach()
  install(FILES pydlr/pydlr.py DESTINATION pydlr)

  add_test(NAME python-test_kernel
    COMMAND ${Python3_EXECUTABLE} -m nose ./test/test_kernel.py -v
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pydlr
    )

  add_test(NAME python-test_convolution
    COMMAND ${Python3_EXECUTABLE} -m nose ./test/test_convolution.py -v
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pydlr
    )

endif()
