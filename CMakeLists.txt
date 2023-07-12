# S2HAT

# Top-level CMakeLists file
# - contains the call to the project() command and some global settings
# - search for external dependencies: MPI, OpenMP, FFTW3
# - add relevant subdirectories to the build

# ----------------------------------------------------------------------------------------#
# General configuration
# ----------------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.22...3.26)
enable_language( Fortran )

# set(CMAKE_OSX_DEPLOYMENT_TARGET "13.1" CACHE STRING "Minimum OS X deployment version")

# Project name and a few useful settings
project(S2HAT
        VERSION 2.56
        DESCRIPTION "Spherical harmonic librairy for high performance CMB data analysis"
        LANGUAGES C Fortran)

# Enforce language standards
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(CMAKE_Fortran_STANDARD 77)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

# We are building libraries that will eventually be linked into shared modules.
# All code should be built with PIC.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# ----- RPATH -----

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if ("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif ("${isSystemDir}" STREQUAL "-1")

# Auxiliary files
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# ----------------------------------------------------------------------------------------#
# External packages
# ----------------------------------------------------------------------------------------#

# Search for MPI 
find_package(MPI REQUIRED)
find_package(MPI COMPONENTS Fortran REQUIRED)

# Search for FFTW
# find_package(FFTW REQUIRED)
find_package(FFTW REQUIRED DOUBLE_LIB)

# Search for HEALPIX
find_package(HEALPIX)

# ----------------------------------------------------------------------------------------#
# Subdirectories
# ----------------------------------------------------------------------------------------#


add_subdirectory(src)

option(BUILD_TESTS "Build the tests" OFF)

if (BUILD_TESTS)
    add_subdirectory(examples)
endif ()
