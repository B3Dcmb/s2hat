# Executable for running mappraiser without TOAST workflow

add_executable(test_run mini_test.c)

set_property(TARGET test_run PROPERTY LINKER_LANGUAGE Fortran)

target_sources(test_run PRIVATE "mini_test.c")
target_include_directories(test_run PRIVATE
        "../include"
        "../lib"
        )


add_compile_options(${MPI_Fortran_COMPILE_FLAGS})

# MPI compile flag for midapack algebra modules
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DMAXCHK=1500 -O3 -fopenmp -fPIC -fno-second-underscore -ffixed-line-length-none -ffree-line-length-none -fallow-argument-mismatch")

# target_compile_definitions(test_run PRIVATE -fallow-argument-mismatch -fopenmp)
target_compile_definitions(test_run PRIVATE -Dx86_64)


option(BUILD_with_FFTW "Build with FFTW" ON)
option(BUILD_with_HEALPIX "Build with HEALPIX" OFF)


if (FFTW_FOUND AND BUILD_with_FFTW)
    #target_compile_definitions(s2hat_std PRIVATE HAVE_FFTW=1)
    target_include_directories(test_run PUBLIC "${FFTW_INCLUDE_DIRS}")
    target_link_libraries(test_run PUBLIC "${FFTW_LIBRARIES}")
    #    if(FFTW_DOUBLE_THREADS_LIB_FOUND)
    #        target_compile_definitions(s2hat_std PRIVATE HAVE_FFTW_THREADS=1)
    #    endif(FFTW_DOUBLE_THREADS_LIB_FOUND)
    target_compile_options(test_run PUBLIC -DFFTW3_C2R)
endif (FFTW_FOUND AND BUILD_with_FFTW)

if (HEALPIX_FOUND AND BUILD_with_HEALPIX)
    #target_compile_definitions(s2hat_std PRIVATE HAVE_HEALPIX=1)
    target_include_directories(test_run PUBLIC "${HEALPIX_INCLUDE_DIRS}")
    target_link_libraries(test_run PUBLIC "${HEALPIX_LIBRARIES}")
    target_compile_options(test_run PUBLIC -DHEALPIX_fft)
endif (HEALPIX_FOUND AND BUILD_with_HEALPIX)


# libraries
target_link_libraries(test_run PRIVATE s2hat_std MPI::MPI_C)
target_link_libraries(test_run PRIVATE MPI::MPI_Fortran)


install(TARGETS test_run 
        DESTINATION ${PROJECT_SOURCE_DIR}/examples/)
