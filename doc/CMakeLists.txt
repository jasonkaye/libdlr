
find_package(Sphinx)

SET(sphinx_top ${CMAKE_CURRENT_BINARY_DIR}/html/index.html)

FILE(GLOB_RECURSE SOURCES *.rst)
EXECUTE_PROCESS(COMMAND sh -c "mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/html")
set(DOC_SOURCE ${CMAKE_CURRENT_BINARY_DIR})

# Create a temporary copy of the doc that will be completed by the automatically generated one.
FILE(GLOB_RECURSE all_rst_files RELATIVE ${CMAKE_SOURCE_DIR}/doc *.rst _* *.py *.png *.cpp Doxyfile)
SET(rst_copy_tar  ${CMAKE_CURRENT_BINARY_DIR}/rst_copy.tar)
add_custom_command (OUTPUT ${rst_copy_tar} DEPENDS ${all_rst_files}
 COMMAND cd ${CMAKE_SOURCE_DIR}/doc && tar cf ${rst_copy_tar} ${all_rst_files} && cd ${CMAKE_BINARY_DIR}/doc && tar xf ${rst_copy_tar} && cp conf.py.patch conf.py )
add_custom_target(rst_copy DEPENDS ${rst_copy_tar})

#add_subdirectory(_doxygen)
add_subdirectory(_doxygen_fortran)

SET(sphinx_cfg ${CMAKE_CURRENT_BINARY_DIR}/conf.py.patch)
file(READ conf.py SPHINXCFG_PATCH)
#string(REPLACE "../c++" "../../c++" SPHINXCFG_PATCH ${SPHINXCFG_PATCH})
string(REPLACE "../src" "../../src" SPHINXCFG_PATCH ${SPHINXCFG_PATCH})
string(REPLACE "../pydlr" "../../pydlr" SPHINXCFG_PATCH ${SPHINXCFG_PATCH})
file(WRITE ${sphinx_cfg} ${SPHINXCFG_PATCH})

add_custom_command(OUTPUT ${sphinx_top} DEPENDS ${SOURCES}
                   COMMAND  export PYTHONPATH=$$PYTHONPATH:${CMAKE_SOURCE_DIR} && 
${SPHINX_EXECUTABLE} -c . -b html ${DOC_SOURCE} html)
add_custom_target(docs_sphinx${DOC_EXT} ALL DEPENDS ${sphinx_top} ${DOC_SOURCE})

add_dependencies(docs_sphinx rst_copy) # we must first have copied the sources
add_dependencies(docs_sphinx docs_doxygen_fortran) 
