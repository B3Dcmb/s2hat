# Executable for running mappraiser without TOAST workflow

add_executable(debug_run mini_test.c)

target_sources(debug_run PRIVATE "mini_test.c")
target_include_directories(debug_run PRIVATE
        "../include"
        "../lib"
        )

# MPI compile flag for midapack algebra modules
# target_compile_definitions(debug_run PRIVATE -fallow-argument-mismatch -fopenmp)

# if (FFTW_FOUND)
#     # #target_compile_definitions(s2hat_std PRIVATE HAVE_FFTW=1)
#     # target_include_directories(debug_run PUBLIC "${FFTW_INCLUDE_DIRS}")
#     # target_link_libraries(debug_run PUBLIC "${FFTW_LIBRARIES}")
#     # #    if(FFTW_DOUBLE_THREADS_LIB_FOUND)
#     # #        target_compile_definitions(s2hat_std PRIVATE HAVE_FFTW_THREADS=1)
#     # #    endif(FFTW_DOUBLE_THREADS_LIB_FOUND)
#     target_compile_definitions(debug_run PUBLIC -DFFTW3_C2R)
# endif (FFTW_FOUND)

# if (HEALPIX_FOUND)
#     # #target_compile_definitions(s2hat_std PRIVATE HAVE_HEALPIX=1)
#     # target_include_directories(debug_run PUBLIC "${HEALPIX_INCLUDE_DIRS}")
#     # target_link_libraries(debug_run PUBLIC "${HEALPIX_LIBRARIES}")
#     # target_compile_options(debug_run PUBLIC -DHEALPIX_fft)
# endif (HEALPIX_FOUND)

# libraries
target_compile_definitions(debug_run PUBLIC -DFFTW3_C2R)
target_link_libraries(debug_run PRIVATE s2hat_std MPI::MPI_C)
target_link_libraries(debug_run PRIVATE s2hat_std MPI::MPI_Fortran)
target_compile_options(debug_run PUBLIC -DHEALPIX_fft)


add_executable(mini_test mini_test.c) # hardcode lib to make it simpler