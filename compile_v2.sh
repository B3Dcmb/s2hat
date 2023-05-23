

# set echo 

# export HEALPIXPATH '-I .... -L ....'    # this is needed for the FFT
export FFTPATH="-I/global/homes/m/mag/software/Healpix_3.82/include -L/global/homes/m/mag/software/Healpix_3.82/"    # this is needed for the FFT
# export DFLAGS='-Dx86_64 -DMAXCHK=1500  -DHEALPIX_fft'
export DFLAGS='-Dx86_64 -DMAXCHK=1500  -DFFTW3_C2R'

FFLAGS='-O3 -std=gnu -fno-second-underscore -ffixed-line-length-none -ffree-line-length-none'

GFORTRAN_VERSION=$(gfortran -v 2>&1 | sed -nE 's/.*?gcc version ([0-9]+)\..*/\1/p')
if [ "${GFORTRAN_VERSION}" -ge "10" ]; then
    FFLAGS="-fallow-argument-mismatch ${FFLAGS}"
fi

export FFLAGS="-O3 -std=gnu -fopenmp -fno-second-underscore -ffixed-line-length-none -ffree-line-length-none -fallow-argument-mismatch"
echo $FFLAGS
# export FFLAGS '-O3 -std=gnu -fno-second-underscore -ffixed-line-length-none -ffree-line-length-none'
export CFLAGS='-O3 -g'

ftn $DFLAGS $FFLAGS  -c s2hat_defs.f90
ftn $DFLAGS $FFLAGS $FFTPATH -c s2hat_types_internal.F90
ftn $DFLAGS $FFLAGS  -c s2hat_pixelization.f90
ftn $DFLAGS $FFLAGS $FFTPATH  -c s2hat_toolbox.F90
ftn $DFLAGS $FFLAGS  -c s2hat_alm2map.f90
ftn $DFLAGS $FFLAGS  -c s2hat_map2alm.F90
ftn $DFLAGS $FFLAGS  -c s2hat_c_interface.f90

cc $DFLAGS $CFLAGS -c s2hat_c_wrappers.c

ar rv /global/u2/m/mag/software/s2hat_cmake_cori/lib/libs2hat_std.a s2hat_alm2map.o s2hat_c_interface.o s2hat_c_wrappers.o s2hat_defs.o s2hat_map2alm.o s2hat_pixelization.o s2hat_toolbox.o s2hat_types_internal.o
cp ../src_fortran/*.mod /global/u2/m/mag/software/s2hat_cmake_cori/lib/
