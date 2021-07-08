

SET(doxygen_top ${CMAKE_CURRENT_BINARY_DIR}/xml/index.xml)

SET(doxygen_cfg ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile.patch)
file(READ Doxyfile DOXCFG_INPUT)
string(REPLACE "../../src" "../../../src" DOXCFG_PATCH ${DOXCFG_INPUT})
file(WRITE ${doxygen_cfg} ${DOXCFG_PATCH})

add_custom_command(OUTPUT ${doxygen_top} DEPENDS ${doxygen_cfg} COMMAND doxygen Doxyfile.patch)
add_custom_target(docs_doxygen_fortran${DOC_EXT} ALL DEPENDS ${doxygen_top} ${doxygen_cfg})

add_dependencies(docs_doxygen_fortran rst_copy)
