# LIBs2hat_std

# add_library(s2hat_std SHARED)
# add_library(s2hat_std STATIC 
#             "${CMAKE_CURRENT_LIST_DIR}/c_wrapper/s2hat_c_wrappers.c" 
#             "${CMAKE_CURRENT_LIST_DIR}/../include/s2hat.h" 
#             "${CMAKE_CURRENT_LIST_DIR}/../include/s2hat_types.h")
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/../include")
# set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/../include")
# target_include_directories(s2hat_std PUBLIC Fortran_MODULE_DIRECTORY)


option(BUILD_with_FFTW "Build with FFTW" ON)
option(BUILD_with_HEALPIX "Build with HEALPIX" OFF)


# add_library(s2hat_std STATIC 
#             "${CMAKE_CURRENT_LIST_DIR}/../include/s2hat.h" 
#             "${CMAKE_CURRENT_LIST_DIR}/../include/s2hat_types.h"
#             "${CMAKE_CURRENT_LIST_DIR}/../include/s2hat_f2c.h"
            # "${CMAKE_CURRENT_LIST_DIR}/../include/s2hat_defs.h")

add_library(s2hat_std STATIC)
# set_property(TARGET s2hat_std PROPERTY LINKER_LANGUAGE Fortran)

# ----------------------------------------------------------------------------------------#
# sources files
# ----------------------------------------------------------------------------------------#

add_subdirectory(src_fortran)
add_subdirectory(c_wrapper)

# ----------------------------------------------------------------------------------------#
# header files
# ----------------------------------------------------------------------------------------#

target_include_directories(s2hat_std PUBLIC
                           "../include"
                           )

# TODO build both a shared and static library? allow to choose?
# set(SOURCE_FILES 
#     "${CMAKE_CURRENT_LIST_DIR}/c_wrapper/s2hat_c_wrappers.c"
#     "s2hat_defs.f90"
#     "s2hat_types_internal.F90"
#     "s2hat_pixelization.f90"
#     "s2hat_toolbox.F90"
#     "s2hat_alm2map.f90"
#     "s2hat_map2alm.F90"
#     "s2hat_c_interface.f90")

# set(SOURCE_HEADERS "${CMAKE_CURRENT_LIST_DIR}/../include/s2hat.h" "${CMAKE_CURRENT_LIST_DIR}/../include/s2hat_types.h")

# target_sources(s2hat_std PUBLIC 
#         "s2hat_c_wrappers.c"
#         "s2hat_defs.f90"
#         "s2hat_types_internal.F90"
#         "s2hat_pixelization.f90"
#         "s2hat_toolbox.F90"
#         "s2hat_alm2map.f90"
#         "s2hat_map2alm.F90"
#         "s2hat_c_interface.f90"
#         "s2hat.h"
#         "s2hat_types.h")


# ----------------------------------------------------------------------------------------#
# compile options
# ----------------------------------------------------------------------------------------#

# target_compile_options(s2hat_std PRIVATE -fno-second-underscore -ffixed-line-length-none -ffree-line-length-none -fallow-argument-mismatch -Wall -Wextra -Wpedantic)
# target_compile_options(s2hat_std PUBLIC -Dx86_64 -DMAXCHK=1500 -O3 -fopenmp -fPIC -fno-second-underscore -ffixed-line-length-none -ffree-line-length-none -fallow-argument-mismatch)
# target_compile_options(s2hat_std PUBLIC -Dx86_64 -DMAXCHK=1500 -O3 -fPIC -fno-second-underscore -ffree-line-length-none -fallow-argument-mismatch)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Dx86_64 -DMAXCHK=1500 -O3 -fopenmp -fPIC -fno-second-underscore -ffixed-line-length-none -ffree-line-length-none -fallow-argument-mismatch")
add_compile_options(${MPI_Fortran_COMPILE_FLAGS})

# target_compile_options(s2hat_std PRIVATE -O3)
            
# ----------------------------------------------------------------------------------------#
# dependencies (searched for in top-level CMakeLists.txt)
# ----------------------------------------------------------------------------------------#

# MPI
# set_target_properties(s2hat_std PROPERTIES COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}")
# set_target_properties(s2hat_std PROPERTIES LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")

target_link_libraries(s2hat_std PUBLIC MPI::MPI_C)
target_link_libraries(s2hat_std PUBLIC MPI::MPI_Fortran)
# target_link_libraries(s2hat_std PUBLIC ${MPI_Fortran_LIBRARIES})


# FFTW

if (FFTW_FOUND AND BUILD_with_FFTW)
    #target_compile_definitions(s2hat_std PRIVATE HAVE_FFTW=1)
    target_include_directories(s2hat_std PUBLIC "${FFTW_INCLUDE_DIRS}")
    target_link_libraries(s2hat_std PUBLIC "${FFTW_LIBRARIES}")
    #    if(FFTW_DOUBLE_THREADS_LIB_FOUND)
    #        target_compile_definitions(s2hat_std PRIVATE HAVE_FFTW_THREADS=1)
    #    endif(FFTW_DOUBLE_THREADS_LIB_FOUND)
    target_compile_options(s2hat_std PUBLIC -DFFTW3_C2R)
endif (FFTW_FOUND AND BUILD_with_FFTW)

if (HEALPIX_FOUND AND BUILD_with_HEALPIX)
    #target_compile_definitions(s2hat_std PRIVATE HAVE_HEALPIX=1)
    target_include_directories(s2hat_std PUBLIC "${HEALPIX_INCLUDE_DIRS}")
    target_link_libraries(s2hat_std PUBLIC "${HEALPIX_LIBRARIES}")
    target_compile_options(s2hat_std PUBLIC -DHEALPIX_fft)
endif (HEALPIX_FOUND AND BUILD_with_HEALPIX)

# if (HEALPIX_FOUND)
# set(HEADER_FILES ${HEALPIXROOT}"/include/chealpix.h")
# target_include_directories(s2hat_std PRIVATE $ENV{HEALPIXROOT}/include)
# target_link_libraries(s2hat_std PRIVATE $ENV{HEALPIXROOT}/lib/libchealpix.a)
# target_compile_definitions(s2hat_std PRIVATE "-DHEALPIXDATA=${HEALPIXROOT}share/healpix/ -ls2hat_std")
# endif ()


# ----------------------------------------------------------------------------------------#
# compile options
# ----------------------------------------------------------------------------------------#

# include(CheckFortranCompilerFlag)
# include(CheckCCompilerFlag)

# macro that checks if flag if supported for C, if so add to s2hat
# macro(add_c_flag_if_avail_s2hat_std FLAG)
#     if (NOT "${FLAG}" STREQUAL "")
#         # create a variable for checking the flag if supported, e.g.:
#         # -fp-model=precise --> c_fp_model_precise
#         string(REGEX REPLACE "^-" "c_" FLAG_NAME "${FLAG}")
#         string(REPLACE "-" "_" FLAG_NAME "${FLAG_NAME}")
#         string(REPLACE " " "_" FLAG_NAME "${FLAG_NAME}")
#         string(REPLACE "=" "_" FLAG_NAME "${FLAG_NAME}")

#         check_c_compiler_flag("${FLAG}" ${FLAG_NAME})

#         if (${FLAG_NAME})
#             target_compile_options(s2hat_std PRIVATE $<$<COMPILE_LANGUAGE:C>:${FLAG}>)
#         endif ()
#     endif ()
# endmacro()

# # macro that checks if list of flags are supported for C
# macro(add_flags_if_avail_s2hat)
#     foreach (FLAG ${ARGN})
#         add_c_flag_if_avail_s2hat("${FLAG}")
#     endforeach ()
# endmacro()

# standard flags for C
# add_flags_if_avail_s2hat("-W" "-Wall" "-Wextra" "-Wshadow")

# ----------------------------------------------------------------------------------------#
# install s2hat_std library in CMAKE_INSTALL_PREFIX
# ----------------------------------------------------------------------------------------#

include(GNUInstallDirs)

set(CMAKE_INSTALL_LIBDIR "${CMAKE_CURRENT_LIST_DIR}/../lib")
# set(CMAKE_INSTALL_FULL_LIBDIR /lib)

# install the library
install(TARGETS s2hat_std
        DESTINATION ${CMAKE_INSTALL_LIBDIR}) # hardcode lib to make it simpler
#         #DESTINATION ${CMAKE_INSTALL_LIBDIR})

# install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY} DESTINATION "${CMAKE_CURRENT_LIST_DIR}/../include")
# install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY} DESTINATION ${CMAKE_INSTALL_LIBDIR})

# install the headers
# install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/
#         DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/s2hat)
