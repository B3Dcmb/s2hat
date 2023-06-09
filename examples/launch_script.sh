S2HATROOT=/global/homes/m/mag/software/s2hat_cmake_fftw/s2hat
# HEALPIXROOT=/global/homes/m/mag/software/Healpix_3.82

# cc s2hat_example.c -fopenmp -O3 -DHEALPIX_fft -g -lgfortran -L$S2HATROOT/lib -I$S2HATROOT/include -ls2hat_std  -o s2hat_example
cc s2hat_example.c -O3 -DFFTW3_C2R -g -lgfortran -L$S2HATROOT/lib -I$S2HATROOT/include -ls2hat_std  -o s2hat_example

cc mini_test.c -lgfortran -L$S2HATROOT/lib -I$S2HATROOT/include -ls2hat_std  -o mini_test
cc examples/mini_test.c -L$S2HATROOT/lib -I$S2HATROOT/include -ls2hat_std  -o examples/mini_test_perl


ftn  s2hat_example.f90 -Wmissing-include-dirs -DFFTW3_C2R  -L$S2HATROOT/lib -I$S2HATROOT/include -ls2hat_std -o s2hat_example


cc mini_test.c -O3 -fopenmp -DFFTW3_C2R -L$S2HATROOT/lib -I$S2HATROOT/include -ls2hat_std  -o mini_test

######



ftn  examples/s2hat_example.f90 -Wmissing-include-dirs -DFFTW3_HC2R  -L$S2HATROOT/lib -I$S2HATROOT/include -ls2hat_std -o examples/s2hat_example_f
cc examples/s2hat_example.c -fopenmp -O3 -DFFTW3_C2R -g -dynamic -lgfortran -L$S2HATROOT/lib -I$S2HATROOT/include -ls2hat_std  -o examples/s2hat_example_c

ftn  examples/s2hat_example.f90 -Wmissing-include-dirs -FFTW3_HC2R  -L$S2HATROOT/lib -I$S2HATROOT/include -ls2hat_std -o examples/s2hat_example_f
cc examples/s2hat_example.c -fopenmp -O3 -DHEALPIX_fft -g -dynamic -lgfortran -L$S2HATROOT/lib -I$S2HATROOT/include -ls2hat_std  -o examples/s2hat_example_c



# S2HATROOT_2=/global/homes/m/mag/software/s2hat_cmake_perlmutter_test/old_install
ftn  s2hat_example.f90 -Wmissing-include-dirs -DFFTW3_HC2R  -L$S2HATROOT_2/lib -I$S2HATROOT_2/src -ls2hat_std -o s2hat_example_f2
cc s2hat_example.c -fopenmp -O3 -DFFTW3_C2R -g -dynamic -lgfortran -L$S2HATROOT_2/lib -I$S2HATROOT_2/src -ls2hat_std  -o s2hat_example_c2
cc mini_test.c -fopenmp -O3 -DFFTW3_C2R -g -dynamic -lgfortran -L$S2HATROOT_2/lib -I$S2HATROOT_2/src -ls2hat_std  -o mini_test_c2

# S2HATROOT_3=/global/homes/m/mag/software/s2hat_cmake_perlmutter_test
ftn  s2hat_example.f90 -Wmissing-include-dirs -DFFTW3_HC2R  -L$S2HATROOT_3/lib -I$S2HATROOT_3/include -ls2hat_std -o s2hat_example_f3
cc s2hat_example.c -fopenmp -O3 -DFFTW3_C2R -g -dynamic -lgfortran -L$S2HATROOT_3/lib -I$S2HATROOT_3/include -ls2hat_std  -o s2hat_example_c3
cc mini_test.c -fopenmp -O3 -DFFTW3_C2R -g -dynamic -lgfortran -L$S2HATROOT_3/lib -I$S2HATROOT_3/include -ls2hat_std  -o mini_test_c3


S2HATROOT_4=/global/homes/m/mag/software/Github_softwares/s2hat
ftn mini_test.c -L$S2HATROOT_4/lib -I$S2HATROOT_4/include -ls2hat_std -fopenmp -O3 -fPIC -Dx86_64 -DMAXCHK=1500 -DFFTW3_C2R -std=gnu99 -MD -MT -g -dynamic  -o mini_test
ftn s2hat_example.c -fopenmp -O3 -DFFTW3_C2R -g -dynamic -L$S2HATROOT_4/lib -I$S2HATROOT_4/include -ls2hat_std  -o s2hat_example_c3
