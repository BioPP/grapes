# CMake script for Grapes
# Author: Julien Dutheil
# Created: 30/05/2018

# Global parameters
cmake_minimum_required(VERSION 3.0.0)
project(grapes CXX)

if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif()

set(CMAKE_CXX_FLAGS "-std=c++11 -O2 -g3 -fmessage-length=0 -Wall")



#static linkage?
if (NOT BUILD_STATIC)
  set(BUILD_STATIC FALSE CACHE BOOL
      "Enable static linkage."
      FORCE)
endif()
if (BUILD_STATIC)
  message(STATUS "Static linkage requested.")
  set (CMAKE_CXX_FLAGS "-static -static-libgcc ${CMAKE_CXX_FLAGS}")
endif()

# Check compression program
# COMPRESS_PROGRAM controls the choice of program
# COMPRESS_EXT can be used to override the file extension
if (NOT COMPRESS_PROGRAM)
  set (COMPRESS_PROGRAM gzip CACHE STRING "Set program for compressing documentation" FORCE)
endif ()
find_program (COMPRESS_BIN NAMES ${COMPRESS_PROGRAM} DOC "${COMPRESS_PROGRAM} compression program")
if (NOT COMPRESS_BIN)
  message (STATUS "${COMPRESS_PROGRAM} program not found, text doc will not be compressed")
else ()
  # Deduce COMPRESS_EXT for known compression programs if not set
  if (NOT COMPRESS_EXT)
    if (${COMPRESS_PROGRAM} STREQUAL "gzip")
      set (COMPRESS_EXT "gz")
    elseif (${COMPRESS_PROGRAM} STREQUAL "bzip2")
      set (COMPRESS_EXT "bz2")
    else ()
      set (COMPRESS_EXT "${COMPRESS_PROGRAM}") # Default: program name (works for xz/lzma)
    endif ()
  endif ()
  # Generate command line args (always add -c to output compressed file to stdout)
  if (${COMPRESS_PROGRAM} STREQUAL "gzip")
    # -n for no timestamp in files (reproducible builds)
    # -9 for maximum compression (lintian error)
    set (COMPRESS_ARGS -c -n -9)
  else ()
    set (COMPRESS_ARGS -c)
  endif ()
  message (STATUS "Found ${COMPRESS_BIN} compression program, using file extension .${COMPRESS_EXT}")
endif ()

# Find dependencies (add install directory to search)
if (CMAKE_INSTALL_PREFIX)
  set (CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}" ${CMAKE_PREFIX_PATH})
endif (CMAKE_INSTALL_PREFIX)

include (GNUInstallDirs)
find_package (bpp-phyl3 1.0.0 REQUIRED)

# Find GSL libraries
find_package( GSL )

if (GSL_FOUND)
  include_directories(${GSL_INCLUDE_DIRS})
  set(LIBS ${LIBS} ${GSL_LIBRARIES})
  message(STATUS "GSL libraries found here:")
  message(STATUS "includes: ${GSL_INCLUDE_DIRS}")
  message(STATUS "dynamic libraries: ${GSL_LIBRARIES}")
endif()

# Set the CMAKE_PREFIX_PATH for the find_library fonction when using non
# standard install location
if (CMAKE_INSTALL_PREFIX)
  set (CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}" ${CMAKE_PREFIX_PATH})
endif (CMAKE_INSTALL_PREFIX)

# Subdirectories
ADD_SUBDIRECTORY(grapes)
#ADD_SUBDIRECTORY(doc)
#ADD_SUBDIRECTORY(man)

# Packager
SET(CPACK_PACKAGE_NAME "grapes")
SET(CPACK_PACKAGE_VENDOR "Julien Y. Dutheil")
SET(CPACK_PACKAGE_VERSION "1.1.0")
SET(CPACK_PACKAGE_VERSION_MAJOR "1")
SET(CPACK_PACKAGE_VERSION_MINOR "1")
SET(CPACK_PACKAGE_VERSION_PATCH "0")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "(Greetings_here_is_yet_another) Rate of Adaptive Protein Evolution Software")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE")
SET(CPACK_RESOURCE_FILE_AUTHORS "${CMAKE_SOURCE_DIR}/AUTHORS")
SET(CPACK_RESOURCE_FILE_INSTALL "${CMAKE_SOURCE_DIR}/INSTALL")
SET(CPACK_SOURCE_GENERATOR "TGZ")
# /!\ This assumes that an external build is used
SET(CPACK_SOURCE_IGNORE_FILES 
       "/build/" 
       "/examples/" 
       "/\\\\.git/" 
       "/\\\\.gitignore" 
       ${CPACK_SOURCE_IGNORE_FILES}
       )

SET(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
SET(CPACK_DEBSOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}_${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}.orig")
INCLUDE(CPack)

#This adds the 'dist' target
ADD_CUSTOM_TARGET(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)

