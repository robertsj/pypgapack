#------------------------------------------------------------------------------#
# pypgapack/CMakeLists.txt
#------------------------------------------------------------------------------#

project(pypgapack)
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
set(DETRAN_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#------------------------------------------------------------------------------#
# OPTIONS
#------------------------------------------------------------------------------#

# Set the default build type (debug)
if( NOT CMAKE_BUILD_TYPE)
  set( CMAKE_BUILD_TYPE Debug CACHE STRING
       "Choose the type of build, options are: Debug Release."
       FORCE)
endif()

# Options for MPI. On or Off. 
SET( PGAPACK_ENABLE_MPI NO CACHE BOOL
  "Enable MPI."
)
# PGAPack optimization.  On or Off.
SET( OPTIMIZE NO CACHE BOOL
  "Enable OPTIMIZE."
)
# PGAPack WL parameter.  32 or 64.
SET( WL 64 CACHE INTEGER
  "Maximum number of fixed source iterations allowed."
)

#------------------------------------------------------------------------------#
# EXTERNAL PACKAGES
#------------------------------------------------------------------------------#

# SWIG for generating Python module
find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH})
set(CMAKE_SWIG_FLAGS "")

# This is a hack for now.  Probably worth generating a CMake finder thing.
IF( PGAPACK_ENABLE_MPI )
  #
ENDIF()

#------------------------------------------------------------------------------#
# INSTALL
#------------------------------------------------------------------------------#

set( CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR} CACHE STRING
     "Installation prefix; default is build directory."
     FORCE
)

#------------------------------------------------------------------------------#
# SOURCE
#------------------------------------------------------------------------------#

include_directories(pgapack/include)
include_directories(pypgapack)
add_subdirectory(pgapack)
add_subdirectory(pypgapack)

#------------------------------------------------------------------------------#
# OUTPUT
#------------------------------------------------------------------------------#

message("++++ System:      ${CMAKE_SYSTEM} ${CMAKE_SYSTEM_PROCESSOR}")
message("++++ C:           ${CMAKE_C_COMPILER}")
message("++++ CXX:         ${CMAKE_CXX_COMPILER}")
message("++++ F90:         ${CMAKE_F90_COMPILER}")
message("++++ Build type:  ${CMAKE_BUILD_TYPE}")
message("++++ Install:     ${CMAKE_INSTALL_PREFIX}")

if( DETRAN_ENABLE_PYTHON )
message("++++ MPI:         enabled")
elseif()
message("++++ MPI:         disabled")
endif()



