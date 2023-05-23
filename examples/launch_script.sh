S2HATROOT=/global/homes/m/mag/software/s2hat_cmake_fftw/s2hat
# HEALPIXROOT=/global/homes/m/mag/software/Healpix_3.82

# cc s2hat_example.c -fopenmp -O3 -DHEALPIX_fft -g -lgfortran -L$S2HATROOT/lib -I$S2HATROOT/include -ls2hat_std  -o s2hat_example
cc s2hat_example.c -O3 -DFFTW3_C2R -g -lgfortran -L$S2HATROOT/lib -I$S2HATROOT/include -ls2hat_std  -o s2hat_example

cc mini_test.c -lgfortran -L$S2HATROOT/lib -I$S2HATROOT/include -ls2hat_std  -o mini_test
cc examples/mini_test.c -L$S2HATROOT/lib -I$S2HATROOT/include -ls2hat_std  -o examples/mini_test_perl


ftn  s2hat_example.f90 -Wmissing-include-dirs -DFFTW3_C2R  -L$S2HATROOT/lib -I$S2HATROOT/include -ls2hat_std -o s2hat_example


cc mini_test.c -O3 -fopenmp -DFFTW3_C2R -L$S2HATROOT/lib -I$S2HATROOT/include -ls2hat_std  -o mini_test
