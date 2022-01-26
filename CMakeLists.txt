# set minimum cmake version
cmake_minimum_required(VERSION 3.18)

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

FILE(GLOB TestPrograms RELATIVE ${CMAKE_SOURCE_DIR} ./test/*.f90)
foreach(f ${TestPrograms})
  get_filename_component(name ${f} NAME)
  add_executable(${name} ${f})
  target_link_libraries(${name} PRIVATE dlr)
  add_test(NAME ${name} COMMAND ${name})
endforeach()

# C-library: dlr_c

file(GLOB dlr_c_SRC CONFIGURE_DEPENDS "./src/dlr_c/*.f*")
add_library(dlr_c SHARED ${dlr_c_SRC})
target_link_libraries(dlr_c PUBLIC dlr BLAS::BLAS LAPACK::LAPACK)
install(TARGETS dlr_c LIBRARY)
set_target_properties(dlr_c PROPERTIES INSTALL_RPATH "@loader_path/")

target_include_directories(dlr_c PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>)
target_include_directories(dlr_c SYSTEM INTERFACE $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/src/dlr_c DESTINATION include FILES_MATCHING PATTERN "*.h")

add_executable(c_example ./test/c_example.c)
target_link_libraries(c_example PRIVATE dlr_c)
add_test(NAME c_example COMMAND c_example)

add_executable(ha_it_c ./test/ha_it.c)
target_link_libraries(ha_it_c PRIVATE dlr_c m)
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
  find_python_module("scipy")

  # Copy pydlr files

  FILE(GLOB_RECURSE pydlr_SRC RELATIVE ${CMAKE_SOURCE_DIR} ./pydlr/*.py)
  foreach(f ${pydlr_SRC})
   configure_file(${f} ${f} COPYONLY)
  endforeach()

  # Install only module files (excluding tests)

  FILE(GLOB pydlr_module_SRC RELATIVE ${CMAKE_SOURCE_DIR} ./pydlr/*.py)
  foreach(f ${pydlr_module_SRC})
   install(FILES ${f} DESTINATION pydlr)
  endforeach()

  # Setup python tests
  
  file(GLOB PythonTests RELATIVE ${CMAKE_SOURCE_DIR} ./pydlr/test/test_*.py)
  foreach(f ${PythonTests})
    get_filename_component(name ${f} NAME)
    add_test(NAME ${name}
      COMMAND ${Python3_EXECUTABLE} -m unittest ${f} -v
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
      )
  endforeach()

endif()

# Documentation

if(BUILD_DOCS)
  add_subdirectory(doc)
endif()
