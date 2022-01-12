# set minimum cmake version
cmake_minimum_required(VERSION 3.12)

# project name and language
project(dlr LANGUAGES Fortran C)

# Path for local cmake modules
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

option(with_python "Install python module: pydlr" OFF)
option(BUILD_SHARED_LIBS "Build shared library, turn OFF for static." ON)
option(POSITION_INDEPENDENT_CODE "Use -fPIC by default" ON)
option(BUILD_DOCS "Build documentation" OFF)

include(CTest)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ../lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ../lib)

# Fortran library: dlr

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY lib)

file(GLOB dlr_SRC CONFIGURE_DEPENDS "./src/*.f*" "./src/utils_external/*.f*" "./src/utils_external/id_dist/*.f*" "./src/utils_external/qlinpack/*.f*")

add_library(dlr ${dlr_SRC})
target_link_libraries(dlr PUBLIC BLAS::BLAS LAPACK::LAPACK)
install(TARGETS dlr LIBRARY)
install(TARGETS dlr ARCHIVE)

# Test programs

add_executable(ha_it ./test/ha_it.f90)
target_link_libraries(ha_it PRIVATE dlr)
add_test(NAME ha_it COMMAND ha_it)

add_executable(ha_mf ./test/ha_mf.f90)
target_link_libraries(ha_mf PRIVATE dlr)
add_test(NAME ha_mf COMMAND ha_mf)

add_executable(ip_exp ./test/ip_exp.f90)
target_link_libraries(ip_exp PRIVATE dlr)
add_test(NAME ip_exp COMMAND ip_exp)

add_executable(conv_exp ./test/conv_exp.f90)
target_link_libraries(conv_exp PRIVATE dlr)
add_test(NAME conv_exp COMMAND conv_exp)

add_executable(dyson_sc ./test/dyson_sc.f90)
target_link_libraries(dyson_sc PRIVATE dlr)
add_test(NAME dyson_sc COMMAND dyson_sc)

# C-library: dlr_c

file(GLOB dlr_c_SRC CONFIGURE_DEPENDS "./src/dlr_c/*.f*")
add_library(dlr_c SHARED ${dlr_c_SRC})
target_link_libraries(dlr_c PUBLIC dlr m BLAS::BLAS LAPACK::LAPACK)
install(TARGETS dlr_c LIBRARY)
set_target_properties(dlr_c PROPERTIES INSTALL_RPATH "@loader_path/")

target_include_directories(dlr_c PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>)
target_include_directories(dlr_c SYSTEM INTERFACE $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/src/dlr_c DESTINATION include FILES_MATCHING PATTERN "*.h")

add_executable(c_example ./test/c_example.c)
target_link_libraries(c_example PRIVATE dlr_c)
add_test(NAME c_example COMMAND c_example)

add_executable(ha_it_c ./test/ha_it.c)
target_link_libraries(ha_it_c PRIVATE dlr_c)
add_test(NAME ha_it_c COMMAND ha_it_c)

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
  install(FILES pydlr/kernel.py DESTINATION pydlr)
  install(FILES pydlr/kernel_fortran.py DESTINATION pydlr)
  install(FILES pydlr/utils.py DESTINATION pydlr)
  install(FILES pydlr/__init__.py DESTINATION pydlr)

  add_test(NAME python-test_kernel
    COMMAND ${Python3_EXECUTABLE} -m nose ./test/test_kernel.py -v
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pydlr
    )

  add_test(NAME python-test_free_greens_function
    COMMAND ${Python3_EXECUTABLE} -m nose ./test/test_free_greens_function.py -v
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pydlr
    )

  add_test(NAME python-test_convolution
    COMMAND ${Python3_EXECUTABLE} -m nose ./test/test_convolution.py -v
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pydlr
    )

  add_test(NAME python-test_dimer
    COMMAND ${Python3_EXECUTABLE} -m nose ./test/test_dimer.py -v
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pydlr
    )

  add_test(NAME python-test_barycentric
    COMMAND ${Python3_EXECUTABLE} -m nose ./test/test_barycentric.py -v
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pydlr
    )

  add_test(NAME python-test_gridding
    COMMAND ${Python3_EXECUTABLE} -m nose ./test/test_gridding.py -v
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pydlr
    )

  add_test(NAME python-test_dyson
    COMMAND ${Python3_EXECUTABLE} -m nose ./test/test_dyson.py -v
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pydlr
    )

  add_test(NAME python-test_semi_circ
    COMMAND ${Python3_EXECUTABLE} -m nose ./test/test_semi_circ.py -v
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pydlr
    )

endif()

# Documentation

if(BUILD_DOCS)
  add_subdirectory(doc)
endif()
