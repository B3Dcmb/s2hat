
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"
#include "sys/types.h"
#include "sys/stat.h"
#include "sys/param.h"
#include "mpi.h"

#include "s2hat.h"

/* this code computes first map2alm transform of some map, collects the alm coefficients in the memory       *
 * of the proc root and outputs them as a binary file 'test_alm.bin'. Then it performs alm2map and collects  *
 * the map in the memory of the proc root. It outputs the map as 'test_map.bin'. Then it computes the power  *
 * spectrum out of the precomputed alm coefficients. It is written as a file 'test_cl.bin'. Note that though *
 * in principle the code performs the 'full loop' map2alm followed by alm2map the outputted map will in not  *
 * agree perfectly with the input due to the bandwidth issues.                                               *
 *                                                                                                           *
 *                                                                                      - rs@apc, 2010/02/10 *
 * Works correctly for a multiple maps and any accepted number of Stokes parameters                          *
 *                                                                                                           *
 *                                                                                      - rs@apc, 2010/09/07 */

int main(int argc, char * argv[])
{
   
   s2hat_int4 pixchoice = PIXCHOICE_HEALPIX;   /* use HEALPIX pixelization */

   s2hat_int4 nside = 64;
   s2hat_int4 nmaps = 1;
   s2hat_int4 lmax = 4*nside-1;
   s2hat_int4 npix=12*nside*nside;
   s2hat_int4 root=0;

   s2hat_int4 myrank, nprocs;

   s2hat_int4 nlmax=lmax, nmmax=lmax;
   s2hat_int4 nstokes=3;  /* no of Stokes params 1, 2 or 3 */
   s2hat_int4 nspec = 1;  /* no of map spectra to be computed - check the rules in the manual */
   s2hat_int4 i, j, lout, mout;
   s2hat_int4 plms=0;     /* i.e., do not assume that plm are precomputed */
   s2hat_int4 nmvals, first_ring, last_ring, map_size;
   s2hat_int4 *mvals;
   s2hat_int8 nplm;
   s2hat_int4 one = 1, zero = 0;

   s2hat_flt8 *local_map, *map;
   s2hat_int4 nrings;
   s2hat_flt8 *local_w8ring;

   s2hat_flt8 *cls;

   FILE *fout;

   s2hat_int4 nalms;
   s2hat_dcomplex *local_alm, *alms;

   s2hat_pixeltype cpixelization;
   s2hat_pixparameters pixpar;
   s2hat_scandef cscan;

   s2hat_flt8 zbounds[2];

//         MPI_Init( &argc, &argv);

//   MPI_Comm_rank( MPI_COMM_WORLD, &myrank);
//   MPI_Comm_size( MPI_COMM_WORLD, &nprocs);

  pixpar.par1 = nside; pixpar.par2 = 0;
  set_pixelization( pixchoice, pixpar, &cpixelization); 

   printf("Hello World \n");
   return( 0);

}
