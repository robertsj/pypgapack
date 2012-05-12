#------------------------------------------------------------------------------#
# pypgapack/pgapack/CMakeLists.txt
#------------------------------------------------------------------------------#

# Configure
message(" WL = ${WL}")
set(PGAPACK_CONFIG_IN  "${CMAKE_CURRENT_SOURCE_DIR}/pgapack_config.h.in")
set(PGAPACK_CONFIG_OUT "${PROJECT_BINARY_DIR}/pgapack/include/pgapack_config.h")
message("PGAPACK_CONFIG_IN = ${PGAPACK_CONFIG_IN}")
configure_file(${PGAPACK_CONFIG_IN} ${PGAPACK_CONFIG_OUT})
include_directories("${PROJECT_BINARY_DIR}/pgapack/include/")

#------------------------------------------------------------------------------#
# SOURCE
#------------------------------------------------------------------------------#

FILE(GLOB PGA_SRC source/*.c)

#------------------------------------------------------------------------------#
# LIBRARY
#------------------------------------------------------------------------------#

ADD_LIBRARY(pgapack
            SHARED
            ${PGA_SRC} 
)

#------------------------------------------------------------------------------#
# LIBRARY
#------------------------------------------------------------------------------#

INSTALL(FILES ${CMAKE_CURRENT_SOURCE_DIR}/include/pgapack.h
        DESTINATION ${CMAKE_INSTALL_PREFIX}/include/pgapack.h)
INSTALL(FILES ${CMAKE_CURRENT_SOURCE_DIR}/include/fakempi.h
        DESTINATION ${CMAKE_INSTALL_PREFIX}/include/fakempi.h)
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/include/pgapack_config.h
        DESTINATION ${CMAKE_INSTALL_PREFIX}/include/pgapack_config.h)
INSTALL_TARGETS(/lib pgapack)
