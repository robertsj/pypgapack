#------------------------------------------------------------------------------#
# pypgapack/pypgapack/CMakeLists.txt
#------------------------------------------------------------------------------#

SET(PYPGAPACK_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#-------------------------------------#
# CONFIGURE
#-------------------------------------#

set(PYPGAPACK_CONFIG_IN  "${CMAKE_CURRENT_SOURCE_DIR}/pypgapack_config.hh.in")
set(PYPGAPACK_CONFIG_OUT "${PROJECT_BINARY_DIR}/pypgapack/pypgapack_config.hh")
message("PYPGAPACK_CONFIG_IN = ${PGAPACK_CONFIG_IN}")
configure_file(${PYPGAPACK_CONFIG_IN} ${PYPGAPACK_CONFIG_OUT})
include_directories("${PROJECT_BINARY_DIR}/pypgapack/")
include_directories("${PROJECT_BINARY_DIR}/pgapack/include/")

#-------------------------------------#
# PYTHON INTERFACE
#-------------------------------------#

SET_SOURCE_FILES_PROPERTIES(pypgapack.i PROPERTIES CPLUSPLUS ON)
SWIG_ADD_MODULE(pypgapack python pypgapack.i )
SWIG_LINK_LIBRARIES(pypgapack
                    pgapack
                    ${PYTHON_LIBRARIES} 
                    ${MPI_LIBRARIES}
                   )
INSTALL(FILES ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
        DESTINATION ${CMAKE_INSTALL_PREFIX}/python/pypgapack)                    
INSTALL(FILES ${CMAKE_BINARY_DIR}/pypgapack/pypgapack.py
        DESTINATION ${CMAKE_INSTALL_PREFIX}/python/pypgapack)
INSTALL(FILES ${CMAKE_BINARY_DIR}/pypgapack/_pypgapack.so
        DESTINATION ${CMAKE_INSTALL_PREFIX}/python/pypgapack)          

