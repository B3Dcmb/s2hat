
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"
#include "sys/types.h"
#include "sys/stat.h"
#include "sys/param.h"
#include "mpi.h"

#include "s2hat.h"

#include "sprng.h"

#define max( x, y) (((x) > (y)) ? (x) : (y))

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
/*                                                                                                           * 
 * formerly s2hat_test.c adapted to test the iterative map2s2hat                                             *
 *                                                                                                           *
 *                                                                                      - rs@cpb, 2022/05/22 */

s2hat_int4 s2hat_map2alm_gausjac( s2hat_int4 plms, s2hat_pixeltype cpixelization, s2hat_scandef cscan, s2hat_int4 nlmax, s2hat_int4 nmmax, s2hat_int4 nmvals, s2hat_int4 *mvals, s2hat_int4 nmaps,
				  s2hat_int4 nstokes, s2hat_int4 first_ring, s2hat_int4 last_ring, s2hat_flt8 *local_w8ring, s2hat_int4 map_size, s2hat_flt8 *local_map, s2hat_int4 lda,
				  s2hat_dcomplex *local_alm, s2hat_int4 nplm, s2hat_flt8 *plm, s2hat_int4 nprocs, s2hat_int4 myrank, s2hat_int4 niter, s2hat_flt8 epsilon,
				  s2hat_int4 *iter_out, s2hat_flt8 *eps_out, MPI_Comm mpi_comm);

s2hat_int4 s2hat_map2alm_cg( s2hat_int4 plms, s2hat_pixeltype cpixelization, s2hat_scandef cscan, s2hat_int4 nlmax, s2hat_int4 nmmax, s2hat_int4 nmvals, s2hat_int4 *mvals, s2hat_int4 nmaps,
			     s2hat_int4 nstokes, s2hat_int4 first_ring, s2hat_int4 last_ring, s2hat_flt8 *local_w8ring, s2hat_int4 map_size, s2hat_flt8 *local_map, s2hat_int4 lda,
			     s2hat_dcomplex *local_alm, s2hat_int4 nplm, s2hat_flt8 *plm, s2hat_int4 nprocs, s2hat_int4 myrank, s2hat_int4 niter, s2hat_flt8 epsilon,
			     s2hat_int4 *iter_out, s2hat_flt8 *eps_out, MPI_Comm mpi_comm);

s2hat_int4 s2hat_map2alm_spin_cg( s2hat_pixeltype cpixelization, s2hat_scandef cscan, s2hat_int4 spin, s2hat_int4 nlmax, s2hat_int4 nmmax, s2hat_int4 nmvals, s2hat_int4 *mvals, s2hat_int4 nmaps,
			          s2hat_int4 first_ring, s2hat_int4 last_ring, s2hat_flt8 *local_w8ring, s2hat_int4 map_size, s2hat_flt8 *local_map, s2hat_int4 lda, s2hat_dcomplex *local_alm,
				  s2hat_int4 nprocs, s2hat_int4 myrank, s2hat_int4 niter, s2hat_flt8 epsilon, s2hat_int4 *iter_out, s2hat_flt8 *eps_out, MPI_Comm mpi_comm);

s2hat_int4 s2hat_map2alm_cg_zero( s2hat_int4 plms, s2hat_pixeltype cpixelization, s2hat_scandef cscan, s2hat_int4 nlmax, s2hat_int4 nmmax, s2hat_int4 nmvals, s2hat_int4 *mvals, s2hat_int4 nmaps,
			          s2hat_int4 nstokes, s2hat_int4 first_ring, s2hat_int4 last_ring, s2hat_flt8 *local_w8ring, s2hat_int4 map_size, s2hat_flt8 *local_map, s2hat_int4 lda,
			          s2hat_dcomplex *local_alm, s2hat_int4 nplm, s2hat_flt8 *plm, s2hat_int4 nprocs, s2hat_int4 myrank, s2hat_int4 niter, s2hat_flt8 epsilon, s2hat_int4 *iter_out,
			          s2hat_flt8 *eps_out, MPI_Comm mpi_comm);

s2hat_flt8 hres_norm( s2hat_int4 nmaps, s2hat_int4 nstokes, s2hat_int4 nlmax, s2hat_int4 nmvals, s2hat_int4 *mvals, s2hat_int4 lda, s2hat_dcomplex *local_alm, MPI_Comm mpi_comm);
s2hat_flt8 mres_norm( s2hat_int4 length, s2hat_flt8 *vect, MPI_Comm mpi_comm);
s2hat_flt8 mres_weighted_norm( s2hat_pixeltype cpixelization, s2hat_scandef cscan, s2hat_int4 nstokes, s2hat_int4 first_ring, s2hat_int4 last_ring, s2hat_flt8 *local_w8ring,
			       s2hat_int4 map_size, s2hat_flt8 *local_map, MPI_Comm mpi_comm);
s2hat_int4 combine_alm_sets( s2hat_int4 nmaps, s2hat_int4 nstokes, s2hat_int4 nlmax, s2hat_int4 nmvals, s2hat_int4 *mvals, s2hat_int4 lda, s2hat_flt8 alpha, s2hat_dcomplex *local_alm1, s2hat_flt8 beta, s2hat_dcomplex *local_alm2, MPI_Comm mpi_comm);

s2hat_int4 get_random_alms( s2hat_int4**, s2hat_int4, s2hat_flt8*, s2hat_int4, s2hat_int4, s2hat_int4*, s2hat_dcomplex*);
void sprng_gaussran( int *stream, int nvect, double *rvect);

s2hat_int4 test = 1;   /* 1 - output aux files; 0 - nah ... */


/* main starts below */

s2hat_int4 main(s2hat_int4 argc, char * argv[])
{

  s2hat_int4 pixchoice = PIXCHOICE_HEALPIX;   /* use HEALPIX pixelization */

  s2hat_int4 nside = 1024;
  s2hat_int4 nmaps = 1;
  s2hat_int4 lmax = (int)(2.6*nside)-1; /*3*nside-1;*/  /*2399;*/ /* 2*nside-1;*/
  s2hat_int4 npix=12*nside*nside;
  s2hat_int4 root=0;

  s2hat_int4 myrank, nprocs;

  s2hat_int4 nlmax=lmax, nlmax1, nmmax=lmax;
  s2hat_int4 nstokes=2;  /* no of Stokes params 1, 2 or 3 */
  s2hat_int4 nspec = 2;  /* no of map spectra to be computed - check the rules in the manual */
  s2hat_int4 i, j, l, m, imap, is, ishft, jshft, lout, mindx, mout, nstreams, spin;
  s2hat_int4 plms=0;     /* i.e., do not assume that plm are precomputed */
  s2hat_int4 nmvals, first_ring, last_ring, map_size;
  s2hat_int4 *mvals;
  s2hat_int8 nplm;
  s2hat_int4 one = 1, zero = 0;

  int iseed = 985456376, sprng_type=SPRNG_LFG;
  
  s2hat_flt8 *local_map, *map;
  s2hat_int4 nrings;
  s2hat_flt8 *local_w8ring;

  s2hat_flt8 *cls;

  char fname[100];
  FILE *fout, *fin;

  s2hat_int4 nalms;
  s2hat_dcomplex *local_alm, *local_alm_inp, *alms;

  s2hat_pixeltype cpixelization;
  s2hat_pixparameters pixpar;
  s2hat_scandef cscan;

  s2hat_flt8 zbounds[2];

  /* the convergence params */
  
  s2hat_int4 *iter_out;
  s2hat_flt8 *alm_res, alm_res_ref, *eps_out, *clamp, *inv_lcut;

  s2hat_int4 niter_max = 1000;
  s2hat_flt8 epsilon = 1.0e-10, done = 1.0, dmone = -1.0, dzero = 0.0;
  
  /* MPI initialization */

  MPI_Init( &argc, &argv);

  MPI_Comm_rank( MPI_COMM_WORLD, &myrank);
  MPI_Comm_size( MPI_COMM_WORLD, &nprocs);

  /* define the maps parameters and data distribution */

  spin =2;
  
  /* - define the pixelization structure -> cpixelization */
  pixpar.par1 = nside; pixpar.par2 = 0;
  set_pixelization( pixchoice, pixpar, &cpixelization); 

  /* - define the scan structure -> cscan : corresponding to the observed sky given by zbounds[0] >= cos theta >= zbounds[1] */

  zbounds[0] = 0.0; zbounds[1] = 0.0;
  zbounds2scan( zbounds, cpixelization, &cscan);

  /* computes the data distribution over the processors */

  get_local_data_sizes( plms, cpixelization, cscan, nlmax, nmmax, myrank, nprocs,
			&nmvals, &first_ring, &last_ring, &map_size, &nplm, root, MPI_COMM_WORLD);  

  /* number of equilatitude rings on this processor */
  nrings = last_ring-first_ring+1;

  /* m-values to be stored on this proc */
  mvals = (int *) calloc( nmvals, sizeof( int));
  find_mvalues( myrank, nprocs, nmmax, nmvals, mvals);

  /* number of alm coefficients stored on this proc */
   nalms = (nlmax+1)*nmvals*nstokes*nmaps;

   /* ----- */
   
  /* generate the alms from the initial power spectrum */

  /* fake initial spectrum for time being */

   nlmax1 = nlmax+1;
   
  cls = (s2hat_flt8 *) calloc( nlmax1*nstokes, sizeof( s2hat_flt8));

  clamp = (s2hat_flt8 *)calloc( nstokes, sizeof( s2hat_flt8));
  inv_lcut = (s2hat_flt8 *)calloc( nstokes, sizeof( s2hat_flt8));  

  clamp[0] = 1.0; inv_lcut[0] = 0.0;
  for(is=1; is<nstokes; is++) { clamp[is] = clamp[is-1]/1.0; inv_lcut[is] += 0.0005; }
  
  for( is=0; is<nstokes; is++)
    for(l=2;l<nlmax1;l++)
      cls[is*nlmax1+l] = clamp[is]*exp(-0.5*(l*inv_lcut[is])*(l*inv_lcut[is]));

  free( clamp); free( inv_lcut);
  
  /*  if( myrank == root) iseed = make_sprng_seed();
      MPI_Bcast( &iseed, 1, MPI_INT, root, MPI_COMM_WORLD);*/
  
  fprintf( stdout, "myrank = %d -> %d\n", myrank, iseed); fflush(stdout);
  
  /* total number of random streams */
  nstreams = 2*nstokes*(nmmax+1);

  int **my_streams;

  my_streams = (int **)calloc( nstokes*nmvals, sizeof( int*));   /* to store pointers to my streams */

  if( zero) {
    /* generate my streams - one per m value */
    for( m=0; m<nmvals; m++) {
       my_streams[m] = init_sprng( sprng_type, mvals[m],  nstreams, iseed, SPRNG_DEFAULT);

       if( mvals[m] == 5) {
          fprintf( stdout, " proc = %d ; m = %d \n", myrank, mvals[m]);
          print_sprng(   my_streams[m]);
          fprintf(stdout,"\n\n");
          for(i=0;i<5;i++) fprintf(stdout, " %.4e,", sprng(my_streams[m]));
          fprintf(stdout, "\n"); fflush(stdout);
       }

       if( mvals[m] == 375) {
          fprintf( stdout, " proc = %d ; m = %d \n", myrank, mvals[m]);
          print_sprng(   my_streams[m]);
          fprintf(stdout,"\n\n");
          for(i=0;i<5;i++) fprintf(stdout, " %.4e,", sprng(my_streams[m]));
          fprintf(stdout, "\n"); fflush(stdout);
       }

     /* fprintf( stdout, "myrank = %d m = %d \n", myrank, mvals[m]); fflush(stdout); */
    }

    fprintf( stdout, " out ! [%d]", myrank); fflush(stdout);
  
    MPI_Barrier( MPI_COMM_WORLD);
  
    MPI_Abort( MPI_COMM_WORLD, 1);

  }

  /* -------- */

  /* define the harmonic coeffcients - the 'ground' truth */


  /* initialize random streams for each m separately */

  fprintf( stdout, "myrank = %d : %d boom\n", myrank, nstreams); fflush(stdout);

  MPI_Barrier( MPI_COMM_WORLD);  
  /*
  int **real_streams, **imag_streams;

  real_streams = (int **)calloc( nstokes*nmvals, sizeof( int*));
  imag_streams = (int **)calloc( nstokes*nmvals, sizeof( int*));
    
  for( is=0; is< nstokes; is++)
    ishft = is*nmvals;
    jshft = is*(nmmax+1);
    for( mindx=0; mindx<nmvals; mindx++) {
      real_streams[ishft+mindx] = init_sprng( sprng_type, 2*(jshft+mvals[mindx]),  nstreams, iseed, SPRNG_DEFAULT);
      imag_streams[ishft+mindx] = init_sprng( sprng_type, 2*(jshft+mvals[mindx])+1,  nstreams, iseed, SPRNG_DEFAULT);    

      * print_sprng(   real_streams[ishft+mindx]); *
      fprintf( stdout, " RAND [%d : %d ] -> \n", is, mindx); * , sprng(real_streams[ishft+mindx]), sprng(imag_streams[ishft+mindx]));* fflush( stdout);

    }
  */
  /* use the streams to generate alms */

  fprintf( stdout, "myrank = %d SPRNG initiated \n", myrank); fflush(stdout);

  MPI_Barrier( MPI_COMM_WORLD);
  
  local_alm_inp = (s2hat_dcomplex *) calloc( nalms, sizeof( s2hat_dcomplex));

  double *revect, *imvect;
  int *rstream;
  
  for( is=0; is< nstokes; is++) {

    fprintf( stdout, "myrank = %d loop Stokes %d \n", myrank, is); fflush(stdout);
 
    jshft = is*(nmmax+1);
    for( mindx=0; mindx<nmvals; mindx++) {

       rstream = init_sprng( sprng_type, 2*(jshft+mvals[mindx]),  nstreams, iseed, SPRNG_DEFAULT);
  
       revect = (double *)calloc( nlmax1-mvals[mindx]+1, sizeof( double));
       /* sprng_gaussran( real_streams[is*nmvals+mindx], nlmax1-mvals[mindx]+1, revect); */
       sprng_gaussran( rstream, nlmax1-mvals[mindx]+1, revect);       

       rstream = init_sprng( sprng_type, 2*(jshft+mvals[mindx])+1,  nstreams, iseed, SPRNG_DEFAULT);       
       imvect = (double *)calloc( nlmax1-mvals[mindx]+1, sizeof( double));
       /* sprng_gaussran( imag_streams[is*nmvals+mindx], nlmax1-mvals[mindx]+1, imvect); */
       sprng_gaussran( rstream, nlmax1-mvals[mindx]+1, imvect);       

       ishft=nlmax1*(is*nmvals+mindx);

       if( mvals[mindx] != 0) {

          for( l=max( mvals[mindx], 2), i=0; l<nlmax1; l++,i++) {
               local_alm_inp[ ishft+l].re = sqrt(0.5*cls[is*nlmax1+l])*revect[i];
               local_alm_inp[ ishft+l].im = sqrt(0.5*cls[is*nlmax1+l])*imvect[i];
          }
	
        } else {
  
          for( l=2, i=0; l<nlmax1; l++, i++) {
             local_alm_inp[ ishft+l].re = sqrt(cls[is*nlmax1+l])*revect[i];
  	     local_alm_inp[ ishft+l].im = 0.0;
          }
        }
       
       free(revect); free(imvect);

       fprintf( stdout, "myrank = %d out %d Stokes %d \n", myrank, mindx, is); fflush(stdout);
 
    }
  }

  free( cls);

  /* free( imag_streams); free(real_streams); */
  
  fprintf( stdout, "alm generation : %d \n", myrank); fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);
  
  cls = (s2hat_flt8 *)calloc( (nlmax+1)*nspec, sizeof( s2hat_flt8));
  collect_cls( nmaps, zero, nstokes, nlmax, nmvals, mvals, nlmax, local_alm_inp, 
               nspec, cls, myrank, nprocs, root, MPI_COMM_WORLD);

  if( myrank == root) 
  {
    sprintf( fname, "test_cl%d_init.bin", spin);
     fout = fopen( fname, "w"); fwrite( cls, sizeof( double), (nlmax+1)*nspec, fout); fclose( fout);
  }

  free( cls);

  MPI_Barrier( MPI_COMM_WORLD);
  
  /* go to the map domain - the standard way */

  /* memory to store the map */
  local_map = (s2hat_flt8 *)calloc( nmaps*map_size*nstokes, sizeof( s2hat_flt8));
   
  s2hat_alm2map_spin( cpixelization, cscan, spin, nlmax, nmmax, nmvals, mvals, nmaps, first_ring, last_ring,
		      map_size, local_map, nlmax, local_alm_inp, nprocs, myrank, MPI_COMM_WORLD);
  
  /* gather complete map in the memory of the proc root */

  if( myrank == root)
  {
      map = (s2hat_flt8 *)calloc( nstokes*npix, sizeof( s2hat_flt8));

      sprintf( fname, "test_map%d.bin", spin);
      fout = fopen( fname,"w"); fclose( fout);
  }

  /* do the loop over maps */
  for( i=0; i<nmaps; i++)
  {
     collect_map(cpixelization,nmaps,i,nstokes,map,first_ring,last_ring,map_size,local_map,myrank,nprocs,root,MPI_COMM_WORLD);

     if( myrank == root)
     {
         fout=fopen( fname,"a"); fwrite( map, sizeof( double), npix*nstokes, fout); fclose( fout);
     }
  }
  
  /* CG iterations */
 
  /* quadrature weights */
  local_w8ring = (s2hat_flt8 *)calloc( nrings*nstokes, sizeof( s2hat_flt8));

  /* i.e. use no weighting */
  for( i=0; i<nrings*nstokes; i++) local_w8ring[i] = 1.0;

  local_alm = (s2hat_dcomplex *) calloc( nalms, sizeof( s2hat_dcomplex));

  /* calculate the alm coefficients - these should match the initial one computed earlier */
  /* - use cojugate gradient solver  */

  iter_out = (s2hat_int4 *)calloc( nmaps, sizeof(s2hat_int4));
  eps_out  = (s2hat_flt8 *)calloc( nmaps, sizeof(s2hat_flt8)); 

  if( (nstokes == 1) || (nstokes == 3)) {
     s2hat_map2alm_cg( plms, cpixelization, cscan, nlmax, nmmax, nmvals, mvals,
  	               nmaps, nstokes, first_ring, last_ring,
	               local_w8ring, map_size, local_map, nlmax,
	               local_alm, nplm, NULL, nprocs, myrank,
		       niter_max, epsilon, iter_out, eps_out,                             /* define the convergence */
		       MPI_COMM_WORLD);
  }

  if( nstokes == 2) {
    s2hat_map2alm_spin_cg( cpixelization, cscan, spin, nlmax, nmmax, nmvals, mvals,
  	                   nmaps, first_ring, last_ring,
	                   local_w8ring, map_size, local_map, nlmax,
	                   local_alm, nprocs, myrank,
		           niter_max, epsilon, iter_out, eps_out,                             /* define the convergence */
		           MPI_COMM_WORLD);
  }

  free( local_w8ring);

  /* output power spectra of the recovered alms */

  cls = (s2hat_flt8 *)calloc( (nlmax+1)*nspec, sizeof( s2hat_flt8));
  collect_cls( nmaps, zero, nstokes, nlmax, nmvals, mvals, nlmax, local_alm, 
               nspec, cls, myrank, nprocs, root, MPI_COMM_WORLD);

  if( myrank == root) 
  {
     sprintf( fname, "test_cl%d_cg.bin", spin);
     fout = fopen( fname, "w"); fwrite( cls, sizeof( double), (nlmax+1)*nspec, fout); fclose( fout);
  }

  free( cls);
  
  /* and compare them with the input ones ... - NB this overwrites the output of the iterations */

  alm_res = (s2hat_flt8 *)calloc( nmaps, sizeof(s2hat_flt8));
  for(imap=0;imap<nmaps;imap++) {
     alm_res_ref=hres_norm( one, nstokes, nlmax, nmvals, mvals, nlmax, &local_alm_inp[imap*nmvals*nlmax1], MPI_COMM_WORLD);    
     combine_alm_sets( one, nstokes, nlmax, nmvals, mvals, nlmax, done, &local_alm[imap*nmvals*nlmax1], -done, &local_alm_inp[imap*nmvals*nlmax1], MPI_COMM_WORLD);
     alm_res[imap]=hres_norm( one, nstokes, nlmax, nmvals, mvals, nlmax, &local_alm[imap*nmvals*nlmax1], MPI_COMM_WORLD);

     alm_res[imap]/=alm_res_ref;     
  }
  
  free( local_alm);
  
  /* report the convergence details */
  
  if( myrank == root) {
    fprintf( stdout,"\n\n");    
    fprintf( stdout, "Conjugate Gradient: \n\n"); fflush(stdout);
    fprintf( stdout, "the number of iterations per map: "); fflush(stdout);
    for(imap=0; imap<nmaps; imap++) fprintf( stdout, " [%d] -> [%d]; ", imap, iter_out[imap]); fflush(stdout);
    fprintf( stdout,"\n\n");
    fprintf( stdout, "the precision reached for each map: "); fflush(stdout);

    for(imap=0; imap<nmaps; imap++) {
      fprintf( stdout, " convergence residuals [%d] -> [%.4e];\n ", imap, eps_out[imap]); fflush(stdout);
      fprintf( stdout, " alm level residuals [%d] -> [%.4e]; ", imap, alm_res[imap]); fflush(stdout);      
    }

    fprintf( stdout,"\n\n");
  }

  free( alm_res);   
  free( iter_out); free(eps_out);

  free( local_alm_inp);

  /* conclude */

  destroy_scan( cscan);
  destroy_pixelization( cpixelization);

  /* and finish up */

  MPI_Finalize();

  return( 0);

}

s2hat_flt8 hres_norm( s2hat_int4 nmaps, s2hat_int4 nstokes, s2hat_int4 nlmax, s2hat_int4 nmvals, s2hat_int4 *mvals, s2hat_int4 lda, s2hat_dcomplex *local_alm, MPI_Comm mpi_comm)

/* computes a proper norm of a complex vector, cvect, of harmonic coefficients *
 *                                                                             *
 * N.B., the HEALPIX convention not tested yet as of 2022/06/02                *
 *                                                            - rs, 2022/05/25 */
  
{
  s2hat_int4 istokes, imap, i, l, m;
  s2hat_flt8 rnorm=0.0, rnorm_loc=0.0;

  /*  s2hat_int4 myrank;

      MPI_Comm_rank( mpi_comm, &myrank); */

  
  if (lda == nlmax) {   /* s2hat convention */
    
    for(imap=i=0; imap<nmaps; imap++)
       for( istokes=0; istokes<nstokes; istokes++)
          for( m=0; m<nmvals; m++)
	    for( l=2, i+=2; l< nlmax+1; l++) {
	        if(mvals[m] !=0) {
	           rnorm_loc += local_alm[i].re*local_alm[i].re+local_alm[i].im*local_alm[i].im;
	        } else {
	           rnorm_loc += 0.5*local_alm[i].re*local_alm[i].re;
	        }
	        i++;           /* global running index */
	    }	 
  } else {
    
    if( lda == nstokes) { /* healpix convention*/

      for(imap=i=0; imap<nmaps; imap++)
            for( m=0; m<nmvals; m++)
	      for( l=2, i+=2*nstokes; l< nlmax+1; l++) {
	          if(mvals[m] !=0) {
		    for( istokes=0; istokes<nstokes; istokes++) {
  	               rnorm_loc += local_alm[i].re*local_alm[i].re+local_alm[i].im*local_alm[i].im;
		       i++;
		    }
	          } else {
		     for( istokes=0; istokes<nstokes; istokes++) {
	                rnorm_loc += 0.5*local_alm[i].re*local_alm[i].re;
		        i++;
		     }
		  }
	       }	 
    } else {
      fprintf( stderr, "wrong lda value, hres_norm"); fflush(stderr);
      MPI_Abort( mpi_comm, 1);
   }
  }

  /*  fprintf( stdout, "[%d] -> %.4e \n", myrank, rnorm_loc); fflush( stdout); */
  
  rnorm_loc *= 2.0; /* include the negative m's as well, NB. m=0 rescaled above */
  
  MPI_Allreduce(&rnorm_loc, &rnorm, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

  /* fprintf( stdout, "[%d] -> %.4e : %.4e \n", myrank, rnorm_loc, rnorm); fflush( stdout); */
  
  rnorm = sqrt( rnorm);
  
  return( rnorm);
}

s2hat_flt8 mres_norm( s2hat_int4 length, s2hat_flt8 *vect, MPI_Comm mpi_comm)

/* computes a norm of a complex vector, cvect, of length elements *
 *                                               - rs, 2022/05/23 */
  
{
  s2hat_int4 i;
  s2hat_flt8 rnorm=0.0, rnorm_loc=0.0;

  for( i=0; i<length; i++) rnorm_loc += vect[i]*vect[i];

  MPI_Allreduce(&rnorm_loc, &rnorm, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

  rnorm = sqrt( rnorm);
  
  return( rnorm);
}

s2hat_flt8 mres_weighted_norm( s2hat_pixeltype cpixelization, s2hat_scandef cscan, s2hat_int4 nstokes, s2hat_int4 first_ring, s2hat_int4 last_ring, s2hat_flt8 *local_w8ring,
			       s2hat_int4 map_size, s2hat_flt8 *local_map, MPI_Comm mpi_comm)

/* computes a weighted dot product of the distributed map         *
 * Tested only for HEALPIX                                        *
 *                                               - rs, 2022/06/02 */
  
{
  s2hat_int4 i, j, ipix, jpix, iring, is;
  s2hat_flt8 rnorm=0.0, rnorm_loc=0.0, rnorth, rsouth, tot_wght;

  for( is = 0; is < nstokes; is++) {
    for( i = is*map_size, j = (is+1)*map_size-1, iring = first_ring; iring < last_ring+1; iring++) {

      tot_wght = cpixelization.parea[iring]*cpixelization.qwght[iring]*local_w8ring[iring-first_ring];

      if(cscan.nfl[iring]) {
         for( rnorth=0.0, ipix=0; ipix < cpixelization.nph[iring]; ipix++, i++) rnorth += local_map[i]*local_map[i];    /* north */
      } else
        i += cpixelization.nph[iring];
    
      if(cscan.sfl[iring]) {
         for( rsouth=0.0, jpix=0; jpix < cpixelization.nph[iring]; jpix++, j--) rsouth += local_map[j]*local_map[j];    /* south */
      } else
        j -= cpixelization.nph[iring];

      rnorm_loc += (rnorth+rsouth)*tot_wght;
    }
  }

  MPI_Allreduce(&rnorm_loc, &rnorm, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

  rnorm = sqrt( rnorm);
  
  return( rnorm);
}

s2hat_int4 combine_alm_sets( s2hat_int4 nmaps, s2hat_int4 nstokes, s2hat_int4 nlmax, s2hat_int4 nmvals, s2hat_int4 *mvals, s2hat_int4 lda, s2hat_flt8 alpha, s2hat_dcomplex *local_alm1, s2hat_flt8 beta, s2hat_dcomplex *local_alm2, MPI_Comm mpi_comm)

/* computes alpha *local_alm1 + beta * local_alm2 and stores it in local_alm1  *
 * local_alm1 and local_alm2 have to be distributed in the same way.           *
 *                                                                             *
 * N.B., the HEALPIX convention not tested yet as of 2022/06/02                *
 *                                                            - rs, 2022/05/31 */
  
{
  s2hat_int4 istokes, imap, i, l, m;
  
  if (lda == nlmax) {   /* s2hat convention */
    
    for(imap=i=0; imap<nmaps; imap++)
       for( istokes=0; istokes<nstokes; istokes++)
          for( m=0; m<nmvals; m++)
	     for( l=0; l< nlmax+1; l++) {
	       
		local_alm1[i].re = alpha*local_alm1[i].re+beta*local_alm2[i].re;
		local_alm1[i].im = alpha*local_alm1[i].im+beta*local_alm2[i].im;		  

	        i++;           /* global running index */
	     }	 
  } else {
    
    if( lda == nstokes) { /* healpix convention*/

      for(imap=i=0; imap<nmaps; imap++)
            for( m=0; m<nmvals; m++)
               for( l=0; l< nlmax+1; l++)
		 for( istokes=0; istokes<nstokes; istokes++) {
		   
		     local_alm1[i].re = alpha*local_alm1[i].re+beta*local_alm2[i].re;
		     local_alm1[i].im = alpha*local_alm1[i].im+beta*local_alm2[i].im;		  

		     i++;           /* global running index */
		 }	 
    } else {
      
      fprintf( stderr, "wrong lda value, hres_norm"); fflush(stderr);
      MPI_Abort( mpi_comm, 1);
      
    }
  }

  return( 1);

}

s2hat_int4 s2hat_map2alm_gausjac( s2hat_int4 plms, s2hat_pixeltype cpixelization, s2hat_scandef cscan, s2hat_int4 nlmax, s2hat_int4 nmmax, s2hat_int4 nmvals, s2hat_int4 *mvals, s2hat_int4 nmaps,
				  s2hat_int4 nstokes, s2hat_int4 first_ring, s2hat_int4 last_ring, s2hat_flt8 *local_w8ring, s2hat_int4 map_size, s2hat_flt8 *local_map, s2hat_int4 lda,
				  s2hat_dcomplex *local_alm, s2hat_int4 nplm, s2hat_flt8 *plm, s2hat_int4 nprocs, s2hat_int4 myrank, s2hat_int4 niter, s2hat_flt8 epsilon, s2hat_int4 *iter_out,
				  s2hat_flt8 *eps_out, MPI_Comm mpi_comm)

/* it is supposed to produce a better estimate of alm coefficients given a pixelized sky map provided in the distributed form. The routine uses fixed point, Gauss-Jacobi solver     *
 * applied for each of the nmaps map separately but for all Stokes parameters. The convergence is reached when the requested number of iterations (niter) or the requested precision * 
 * (epsilon) is reached. The actually reached values of these parameters are stored in iter_out and eps_out for each processed map.                                                  *
 * This version starts with the initial guess set to the pseudo-alms. The computational cost is one map2alm transform to compute the rhs and the initial guess and then one alm2map  *
 * followed by one map2alm per iteration.                                                                                                                                            *
 *                                                                                                                                                                                   *
 * The total cost is then :        (niter+1) of alm2map + niter of map2alm                                                                                       - rs@cpb 2022/05/23 */ 

{
  s2hat_int4 one = 1, root = 0, zero = 0;
  s2hat_int4 i, imap, ishft, iter, nalms, nspec;
  s2hat_flt8 done = 1.0, dmone = -1.0, eps, eps_ref;
  s2hat_flt8 *cls, *local_map0;
  s2hat_dcomplex *local_alm0;

  char fname[50];
  FILE *fout;
  
  /* 0th order, non-iterative solution */
  
  s2hat_map2alm( plms, cpixelization, cscan, nlmax, nmmax, nmvals, mvals, nmaps, nstokes, first_ring, last_ring,
	         local_w8ring, map_size, local_map, nlmax, local_alm, nplm, plm, nprocs, myrank, mpi_comm);

  /* if iterations requested */

  if( niter != 0) {
    
    nspec=nstokes; /* compute only auto-spectra */
    
    for(imap=0; imap< nmaps; imap++) {    /* for each map separately */

      /* eps_ref=hres_norm( one, nstokes, nlmax, nmvals, mvals, nlmax, &local_alm[imap*nalms], mpi_comm); */

       eps_ref = mres_norm( map_size*nstokes, &local_map[imap*nstokes*map_size], mpi_comm);
       
       eps = 10.0*epsilon*eps_ref;
       nalms = nstokes*nmvals*(nlmax+1);  /* number of alm coaffs on this proc per single maps */

       local_map0 = (s2hat_flt8 *)calloc( map_size*nstokes, sizeof( s2hat_flt8));
       local_alm0 = (s2hat_dcomplex *)calloc( nalms, sizeof( s2hat_dcomplex));

       for( iter=0; (iter < niter) && (eps > epsilon*eps_ref); iter++) {

           s2hat_alm2map( plms, cpixelization, cscan, nlmax, nmmax, nmvals, mvals, one, nstokes, first_ring, last_ring,
   	  	          map_size, local_map0, nlmax, &local_alm[imap*nalms], nplm, plm, nprocs, myrank, mpi_comm);
	   
	   ishft = imap*map_size*nstokes;
	   for(i=0; i<map_size*nstokes; i++) local_map0[i] -= local_map[ishft+i];    /* just a single, selected map */

           /* compute map-level residual */

	   eps = mres_norm( map_size*nstokes, local_map0, mpi_comm);
	   
           s2hat_map2alm( plms, cpixelization, cscan, nlmax, nmmax, nmvals, mvals, one, nstokes, first_ring, last_ring,
               	          local_w8ring, map_size, local_map0, nlmax, local_alm0, nplm, plm, nprocs, myrank, mpi_comm);

           /* update alms */

	   ishft = imap*nalms;
  	   /* for( i=0; i< nalms; i++) { local_alm[ishft+i].re -= local_alm0[i].re; local_alm[ishft+i].im -= local_alm0[i].im; }     this should be coded more precisely */

           combine_alm_sets( one, nstokes, nlmax, nmvals, mvals, nlmax, done, &local_alm[ishft], dmone, local_alm0, mpi_comm);
	   
           if( test) {

              cls = (s2hat_flt8 *)calloc( (nlmax+1)*nspec, sizeof( s2hat_flt8));
              collect_cls( nmaps, imap, nstokes, nlmax, nmvals, mvals, nlmax, local_alm, 
                           nspec, cls, myrank, nprocs, root, mpi_comm);

	      if( myrank == root) {
	         sprintf(fname, "gj_cls_iter_%d_map_%d.bin", iter, imap);
	         fout = fopen( fname, "w"); fwrite( cls, sizeof( s2hat_flt8), (nlmax+1)*nspec, fout); fclose( fout);
	      }
	      
	      free(cls);
	   }
	   
           /* compute the increment size the same on all procs */

           /* eps=hres_norm( one, nstokes, nlmax, nmvals, mvals, nlmax, local_alm0, mpi_comm); */

	   if( myrank == root) { fprintf( stdout, "GJ: finishing iteration : %d - map [%d], precision %.4e \n", iter, imap, eps/eps_ref); fflush( stdout); }

       }  /* over iterations */

       free( local_map0);
       free( local_alm0);
  
       iter_out[imap] = iter;
       eps_out[imap] = eps/eps_ref;
    }
  } /* if iterations are performed */
     
  return( 1);

}

s2hat_int4 s2hat_map2alm_cg( s2hat_int4 plms, s2hat_pixeltype cpixelization, s2hat_scandef cscan, s2hat_int4 nlmax, s2hat_int4 nmmax, s2hat_int4 nmvals, s2hat_int4 *mvals, s2hat_int4 nmaps,
			     s2hat_int4 nstokes, s2hat_int4 first_ring, s2hat_int4 last_ring, s2hat_flt8 *local_w8ring, s2hat_int4 map_size, s2hat_flt8 *local_map, s2hat_int4 lda,
			     s2hat_dcomplex *local_alm, s2hat_int4 nplm, s2hat_flt8 *plm, s2hat_int4 nprocs, s2hat_int4 myrank, s2hat_int4 niter, s2hat_flt8 epsilon, s2hat_int4 *iter_out,
			     s2hat_flt8 *eps_out, MPI_Comm mpi_comm)

/* it is supposed to produce a better estimate of alm coefficients given a pixelized sky map provided in the distributed form. The routine uses Conjugate Gradient solver            *
 * applied for each of the nmaps map separately but for all Stokes parameters. The convergence is reached when the requested number of iterations (niter) or the requested precision * 
 * (epsilon) is reached. The actually reached values of these parameters are stored in iter_out and eps_out for each processed map.                                                  *
 * This version starts with the initial guess set to pseudo-alms. The computational cost is one map2alm transform to compute the rhs and one alm2map followed by map2alm to compute  *
 * initial guess and of then one alm2map followed by one map2alm per iteration.                                                                                                      *
 *                                                                                                                                                                                   *
 * The total cost is then :        (niter+2) of alm2map + (niter+1) of map2alm                                                                                   - rs@cpb 2022/05/23 */ 

{

  s2hat_int4 one = 1, zero = 0, root = 0;
  s2hat_int4 i, imap, ishft, jshft, iter, nalms, nspec;
  s2hat_flt8 done = 1.0, dzero=0.0, eps, epsold, eps_ref, alpha_k, beta_k, ypk_norm, tmp;
  s2hat_flt8 *cls, *local_map0;
  s2hat_dcomplex *local_alm0, *local_alm1, *local_alm2;

  char fname[50];
  FILE *fout;
  
  /* 0th order, non-iterative solution = rhs of the equation, Ax = b and the initial guess, x0 */
  
  s2hat_map2alm( plms, cpixelization, cscan, nlmax, nmmax, nmvals, mvals, nmaps, nstokes, first_ring, last_ring,
	         local_w8ring, map_size, local_map, nlmax, local_alm, nplm, plm, nprocs, myrank, mpi_comm);

  /* if iterations requested */

  if( niter != 0) {

     nspec=nstokes; /* compute only auto-spectra */
    
     nalms = nstokes*nmvals*(nlmax+1); /* number of alm coaffs on this proc */

     for( imap=0; imap<nmaps; imap++) {

       	ishft = imap*nalms;
	jshft = imap*map_size*nstokes;
	
        local_map0 = (s2hat_flt8 *)calloc( map_size*nstokes, sizeof( s2hat_flt8));

        /* compute r0 = p0 */
        s2hat_alm2map( plms, cpixelization, cscan, nlmax, nmmax, nmvals, mvals, one, nstokes, first_ring, last_ring,
  	               map_size, local_map0, nlmax, &local_alm[ishft], nplm, plm, nprocs, myrank, mpi_comm);

        for(i=0; i<map_size*nstokes; i++) local_map0[i] = local_map[jshft+i]-local_map0[i];

        local_alm0 = (s2hat_dcomplex *)calloc( nalms, sizeof( s2hat_dcomplex));    /* stores r_k */
     
        s2hat_map2alm( plms, cpixelization, cscan, nlmax, nmmax, nmvals, mvals, one, nstokes, first_ring, last_ring,
  	               local_w8ring, map_size, local_map0, nlmax, local_alm0, nplm, plm, nprocs, myrank, mpi_comm);

        eps_ref=hres_norm( one, nstokes, nlmax, nmvals, mvals, nlmax, &local_alm[ishft], mpi_comm);

	eps=hres_norm( one, nstokes, nlmax, nmvals, mvals, nlmax, local_alm0, mpi_comm);

	if( myrank == root) { fprintf( stdout, "PCG: initiating iteration : map [%d], precision %.4e \n", imap, eps/eps_ref); fflush( stdout); }
	
        local_alm1 = (s2hat_dcomplex *)calloc( nalms, sizeof( s2hat_dcomplex));    /* stores p_k */

        combine_alm_sets( one, nstokes, nlmax, nmvals, mvals, nlmax, dzero, local_alm1, done, local_alm0, mpi_comm);  /* copy r_0 to p_0 */
	
        local_alm2 = (s2hat_dcomplex *)calloc( nalms, sizeof( s2hat_dcomplex));   /* stores A p_k */
     
        for( iter=0; (iter < niter) && (eps > epsilon*eps_ref); iter++) {

           /* pk -> Y pk */
       
           s2hat_alm2map( plms, cpixelization, cscan, nlmax, nmmax, nmvals, mvals, one, nstokes, first_ring, last_ring,
   	 	          map_size, local_map0, nlmax, local_alm1, nplm, plm, nprocs, myrank, mpi_comm);

           ypk_norm = mres_weighted_norm( cpixelization, cscan, nstokes, first_ring, last_ring, local_w8ring, map_size, local_map0, mpi_comm);

           /* ypk_norm = mres_norm( map_size*nstokes, local_map0, mpi_comm)*sqrt(cpixelization.parea[0]);    * this works for equal-area pixelizations only - e.g., HEALPIX */
	   
	   alpha_k = eps*eps/ypk_norm/ypk_norm;  /* if( myrank == root) { fprintf( stdout, "alpha = %.4e eps = %.4e ypk = %.4e [%.4e] \n", alpha_k, eps, ypk_norm, tmp); fflush(stdout); }*/

	   /* x_k+1 */

           combine_alm_sets( one, nstokes, nlmax, nmvals, mvals, nlmax, done, &local_alm[ishft], alpha_k, local_alm1, mpi_comm);
	   
	   /* r_k+1 */
	 
           s2hat_map2alm( plms, cpixelization, cscan, nlmax, nmmax, nmvals, mvals, one, nstokes, first_ring, last_ring,
  	                 local_w8ring, map_size, local_map0, nlmax, local_alm2, nplm, plm, nprocs, myrank, mpi_comm);

	   /* tmp=hres_norm( one, nstokes, nlmax, nmvals, mvals, nlmax, local_alm0, mpi_comm); if( myrank == root) { fprintf(stdout," |p_k| = %.4e \n", tmp); fflush( stdout); } */
	   
	   /* for( i=0; i< nalms; i++) { local_alm0[i].re -= alpha_k*local_alm2[i].re; local_alm0[i].im -= alpha_k*local_alm2[i].im; } */

           combine_alm_sets( one, nstokes, nlmax, nmvals, mvals, nlmax, done, local_alm0, -alpha_k, local_alm2, mpi_comm);

           if( test) {

              cls = (s2hat_flt8 *)calloc( (nlmax+1)*nspec, sizeof( s2hat_flt8));
              collect_cls( nmaps, imap, nstokes, nlmax, nmvals, mvals, nlmax, local_alm, 
                           nspec, cls, myrank, nprocs, root, mpi_comm);

	      if( myrank == root) {
 	         sprintf(fname, "cg0_cls_iter_%d_map_%d.bin", iter, imap);
	         fout = fopen( fname, "w"); fwrite( cls, sizeof( s2hat_flt8), (nlmax+1)*nspec, fout); fclose( fout);
	      }
	      
	      free(cls);
	   }
	 
         /* compute the error the same on all procs */

	 epsold = eps;
         eps=hres_norm( one, nstokes, nlmax, nmvals, mvals, nlmax, local_alm0, mpi_comm);

         if( myrank == root) { fprintf( stdout, "PCG: finishing iteration : %d - map [%d], precision %.4e \n", iter, imap, eps/eps_ref); fflush( stdout); }
	 
	 if( eps < epsilon*eps_ref) break;
	 
         beta_k = eps*eps/epsold/epsold;
	 
         /* p_k+1 */
       
	 /* for( i=0; i< nalms; i++) { local_alm1[i].re = local_alm0[i].re+beta_k*local_alm1[i].re; local_alm1[i].im = local_alm0[i].im+beta_k*local_alm1[i].im; } */

         combine_alm_sets( one, nstokes, nlmax, nmvals, mvals, nlmax, beta_k, local_alm1, done, local_alm0, mpi_comm);  /* output -> local_alm1 */
	 
       }  /* over iterations */

       free( local_map0);
       free( local_alm0); free( local_alm1); free( local_alm2);
  
       iter_out[imap] = iter;
       eps_out[imap]  = eps/eps_ref;

     }   /* over maps */
     
  } /* if iterations are performed */
     
  return( 1);

}

s2hat_int4 s2hat_map2alm_spin_cg( s2hat_pixeltype cpixelization, s2hat_scandef cscan, s2hat_int4 spin, s2hat_int4 nlmax, s2hat_int4 nmmax, s2hat_int4 nmvals, s2hat_int4 *mvals, s2hat_int4 nmaps,
			          s2hat_int4 first_ring, s2hat_int4 last_ring, s2hat_flt8 *local_w8ring, s2hat_int4 map_size, s2hat_flt8 *local_map, s2hat_int4 lda, s2hat_dcomplex *local_alm,
				  s2hat_int4 nprocs, s2hat_int4 myrank, s2hat_int4 niter, s2hat_flt8 epsilon, s2hat_int4 *iter_out, s2hat_flt8 *eps_out, MPI_Comm mpi_comm)

/* it is supposed to produce a better estimate of alm coefficients given a pixelized sky map provided in the distributed form. The routine uses Conjugate Gradient solver            *
 * applied for each of the nmaps map separately but for all Stokes parameters. The convergence is reached when the requested number of iterations (niter) or the requested precision * 
 * (epsilon) is reached. The actually reached values of these parameters are stored in iter_out and eps_out for each processed map.                                                  *
 * This version starts with the initial guess set to pseudo-alms. The computational cost is one map2alm transform to compute the rhs and one alm2map followed by map2alm to compute  *
 * initial guess and of then one alm2map followed by one map2alm per iteration.                                                                                                      *
 *                                                                                                                                                                                   *
 * The total cost is then :        (niter+2) of alm2map + (niter+1) of map2alm                                                                                   - rs@cpb 2022/05/23 *
 *                                                                                                                                                                                   *
 * Adapted to the spin transforms                                                                                                                                - rs@cpb 2022/06/06 */  

{

  s2hat_int4 one = 1, zero = 0, root = 0;
  s2hat_int4 i, imap, ishft, jshft, iter, nalms, nspec, nstokes=2;
  s2hat_flt8 done = 1.0, dzero=0.0, eps, epsold, eps_ref, alpha_k, beta_k, ypk_norm, tmp;
  s2hat_flt8 *cls, *local_map0;
  s2hat_dcomplex *local_alm0, *local_alm1, *local_alm2;

  char fname[50];
  FILE *fout;
  
  /* 0th order, non-iterative solution = rhs of the equation, Ax = b and the initial guess, x0 */
  
  s2hat_map2alm_spin( cpixelization, cscan, spin, nlmax, nmmax, nmvals, mvals, nmaps,first_ring, last_ring,
	              local_w8ring, map_size, local_map, nlmax, local_alm, nprocs, myrank, mpi_comm);

  nspec=nstokes; /* compute only auto-spectra */
  
  if( test) {

    for( imap=0; imap<nmaps; imap++) {
       cls = (s2hat_flt8 *)calloc( (nlmax+1)*nspec, sizeof( s2hat_flt8));
       collect_cls( nmaps, imap, nstokes, nlmax, nmvals, mvals, nlmax, local_alm, 
                    nspec, cls, myrank, nprocs, root, mpi_comm);

       if( myrank == root) {
          sprintf(fname, "cg0_cls%d_niter_map_%d.bin", spin, imap);
          fout = fopen( fname, "w"); fwrite( cls, sizeof( s2hat_flt8), (nlmax+1)*nspec, fout); fclose( fout);
       }
       free( cls);
    }
  }
	      
  /* if iterations requested */

  if( niter != 0) {
    
     nalms = nstokes*nmvals*(nlmax+1); /* number of alm coaffs on this proc */

     for( imap=0; imap<nmaps; imap++) {

       	ishft = imap*nalms;
	jshft = imap*map_size*nstokes;
	
        local_map0 = (s2hat_flt8 *)calloc( map_size*nstokes, sizeof( s2hat_flt8));

        /* compute r0 = p0 */
        s2hat_alm2map_spin( cpixelization, cscan, spin, nlmax, nmmax, nmvals, mvals, one, first_ring, last_ring,
  	               map_size, local_map0, nlmax, &local_alm[ishft], nprocs, myrank, mpi_comm);

        for(i=0; i<map_size*nstokes; i++) local_map0[i] = local_map[jshft+i]-local_map0[i];

        local_alm0 = (s2hat_dcomplex *)calloc( nalms, sizeof( s2hat_dcomplex));    /* stores r_k */
     
        s2hat_map2alm_spin( cpixelization, cscan, spin, nlmax, nmmax, nmvals, mvals, one, first_ring, last_ring,
  	               local_w8ring, map_size, local_map0, nlmax, local_alm0, nprocs, myrank, mpi_comm);

        eps_ref=hres_norm( one, nstokes, nlmax, nmvals, mvals, nlmax, &local_alm[ishft], mpi_comm);

	eps=hres_norm( one, nstokes, nlmax, nmvals, mvals, nlmax, local_alm0, mpi_comm);

	if( myrank == root) { fprintf( stdout, "PCG: initiating iteration : map [%d], precision %.4e \n", imap, eps/eps_ref); fflush( stdout); }
	
        local_alm1 = (s2hat_dcomplex *)calloc( nalms, sizeof( s2hat_dcomplex));    /* stores p_k */

        combine_alm_sets( one, nstokes, nlmax, nmvals, mvals, nlmax, dzero, local_alm1, done, local_alm0, mpi_comm);  /* copy r_0 to p_0 */
	
        local_alm2 = (s2hat_dcomplex *)calloc( nalms, sizeof( s2hat_dcomplex));   /* stores A p_k */
     
        for( iter=0; (iter < niter) && (eps > epsilon*eps_ref); iter++) {

           /* pk -> Y pk */
       
	  s2hat_alm2map_spin( cpixelization, cscan, spin, nlmax, nmmax, nmvals, mvals, one, first_ring, last_ring,
   	 	          map_size, local_map0, nlmax, local_alm1, nprocs, myrank, mpi_comm);

           ypk_norm = mres_weighted_norm( cpixelization, cscan, nstokes, first_ring, last_ring, local_w8ring, map_size, local_map0, mpi_comm);

           /* ypk_norm = mres_norm( map_size*nstokes, local_map0, mpi_comm)*sqrt(cpixelization.parea[0]);    * this works for equal-area pixelizations only - e.g., HEALPIX */
	   
	   alpha_k = eps*eps/ypk_norm/ypk_norm;  /* if( myrank == root) { fprintf( stdout, "alpha = %.4e eps = %.4e ypk = %.4e [%.4e] \n", alpha_k, eps, ypk_norm, tmp); fflush(stdout); }*/

	   /* x_k+1 */

           combine_alm_sets( one, nstokes, nlmax, nmvals, mvals, nlmax, done, &local_alm[ishft], alpha_k, local_alm1, mpi_comm);
	   
	   /* r_k+1 */
	 
           s2hat_map2alm_spin( cpixelization, cscan, spin, nlmax, nmmax, nmvals, mvals, one, first_ring, last_ring,
  	                 local_w8ring, map_size, local_map0, nlmax, local_alm2, nprocs, myrank, mpi_comm);

	   /* tmp=hres_norm( one, nstokes, nlmax, nmvals, mvals, nlmax, local_alm0, mpi_comm); if( myrank == root) { fprintf(stdout," |p_k| = %.4e \n", tmp); fflush( stdout); } */
	   
	   /* for( i=0; i< nalms; i++) { local_alm0[i].re -= alpha_k*local_alm2[i].re; local_alm0[i].im -= alpha_k*local_alm2[i].im; } */

           combine_alm_sets( one, nstokes, nlmax, nmvals, mvals, nlmax, done, local_alm0, -alpha_k, local_alm2, mpi_comm);

           if( test) {

              cls = (s2hat_flt8 *)calloc( (nlmax+1)*nspec, sizeof( s2hat_flt8));
              collect_cls( nmaps, imap, nstokes, nlmax, nmvals, mvals, nlmax, local_alm, 
                           nspec, cls, myrank, nprocs, root, mpi_comm);

	      if( myrank == root) {
		sprintf(fname, "cg0_cls%d_iter_%d_map_%d.bin", spin, iter, imap);
	         fout = fopen( fname, "w"); fwrite( cls, sizeof( s2hat_flt8), (nlmax+1)*nspec, fout); fclose( fout);
	      }
	      
	      free(cls);
	   }
	 
         /* compute the error the same on all procs */

	 epsold = eps;
         eps=hres_norm( one, nstokes, nlmax, nmvals, mvals, nlmax, local_alm0, mpi_comm);

         if( myrank == root) { fprintf( stdout, "PCG: finishing iteration : %d - map [%d], precision %.4e \n", iter, imap, eps/eps_ref); fflush( stdout); }
	 
	 if( eps < epsilon*eps_ref) break;
	 
         beta_k = eps*eps/epsold/epsold;
	 
         /* p_k+1 */
       
	 /* for( i=0; i< nalms; i++) { local_alm1[i].re = local_alm0[i].re+beta_k*local_alm1[i].re; local_alm1[i].im = local_alm0[i].im+beta_k*local_alm1[i].im; } */

         combine_alm_sets( one, nstokes, nlmax, nmvals, mvals, nlmax, beta_k, local_alm1, done, local_alm0, mpi_comm);  /* output -> local_alm1 */
	 
       }  /* over iterations */

       free( local_map0);
       free( local_alm0); free( local_alm1); free( local_alm2);
  
       iter_out[imap] = iter;
       eps_out[imap]  = eps/eps_ref;

     }   /* over maps */
     
  } /* if iterations are performed */
     
  return( 1);

}

s2hat_int4 s2hat_map2alm_cg_zero( s2hat_int4 plms, s2hat_pixeltype cpixelization, s2hat_scandef cscan, s2hat_int4 nlmax, s2hat_int4 nmmax, s2hat_int4 nmvals, s2hat_int4 *mvals, s2hat_int4 nmaps,
			          s2hat_int4 nstokes, s2hat_int4 first_ring, s2hat_int4 last_ring, s2hat_flt8 *local_w8ring, s2hat_int4 map_size, s2hat_flt8 *local_map, s2hat_int4 lda,
			          s2hat_dcomplex *local_alm, s2hat_int4 nplm, s2hat_flt8 *plm, s2hat_int4 nprocs, s2hat_int4 myrank, s2hat_int4 niter, s2hat_flt8 epsilon, s2hat_int4 *iter_out,
			          s2hat_flt8 *eps_out, MPI_Comm mpi_comm)

/* it is supposed to produce a better estimate of alm coefficients given a pixelized sky map provided in the distributed form. The routine uses Conjugate Gradient solver            *
 * applied for each of the nmaps map separately but for all Stokes parameters. The convergence is reached when the requested number of iterations (niter) or the requested precision * 
 * (epsilon) is reached. The actually reached values of these parameters are stored in iter_out and eps_out for each processed map.                                                  *
 * This version starts with the initial guess set to zero. The computational cost is one map2alm transform to compute the rhs of the equation and then one alm2map followed by one   *
 * map2alm per iteration.                                                                                                                                                            *
 *                                                                                                                                                                                   *
 * The total cost is then :        (niter+1) of alm2map + niter of map2alm                                                                                       - rs@cpb 2022/06/03 */ 

{

  s2hat_int4 one = 1, zero = 0, root = 0;
  s2hat_int4 i, imap, ishft, jshft, iter, nalms, nspec;
  s2hat_flt8 done = 1.0, dzero=0.0, eps, epsold, eps_ref, alpha_k, beta_k, ypk_norm, tmp;
  s2hat_flt8 *cls, *local_map0;
  s2hat_dcomplex *local_alm0, *local_alm1, *local_alm2;

  char fname[50];
  FILE *fout;


  /* 0th order, non-iterative solution = rhs of the equation, Ax = b and the initial guess, x0 */
  
  s2hat_map2alm( plms, cpixelization, cscan, nlmax, nmmax, nmvals, mvals, nmaps, nstokes, first_ring, last_ring,
	         local_w8ring, map_size, local_map, nlmax, local_alm, nplm, plm, nprocs, myrank, mpi_comm);

  if( niter != 0) {    /* if iterations requested */

     nspec=nstokes; /* compute only auto-spectra */
    
     nalms = nstokes*nmvals*(nlmax+1); /* number of alm coaffs on this proc */

     for( imap=0; imap<nmaps; imap++) {

       	ishft = imap*nalms;
	jshft = imap*map_size*nstokes;
	
        local_map0 = (s2hat_flt8 *)calloc( map_size*nstokes, sizeof( s2hat_flt8));

        /* compute r0 = p0 = rhs = local_alm*/

        local_alm0 = (s2hat_dcomplex *)calloc( nalms, sizeof( s2hat_dcomplex));    /* stores r_k */
     
        combine_alm_sets( one, nstokes, nlmax, nmvals, mvals, nlmax, dzero, local_alm0, done, &local_alm[ishft], mpi_comm);  /* copy b to r_0 */

        eps = eps_ref = hres_norm( one, nstokes, nlmax, nmvals, mvals, nlmax, local_alm0, mpi_comm);

	if( myrank == root) { fprintf( stdout, "PCG: initiating iteration : map [%d], precision %.4e \n", imap, eps/eps_ref); fflush( stdout); }
	
        local_alm1 = (s2hat_dcomplex *)calloc( nalms, sizeof( s2hat_dcomplex));    /* stores p_k */

        combine_alm_sets( one, nstokes, nlmax, nmvals, mvals, nlmax, dzero, local_alm1, done, local_alm0, mpi_comm);  /* copy r_0 to p_0 */
	
        local_alm2 = (s2hat_dcomplex *)calloc( nalms, sizeof( s2hat_dcomplex));   /* stores A p_k */

        combine_alm_sets( one, nstokes, nlmax, nmvals, mvals, nlmax, dzero, &local_alm[ishft], done, local_alm2, mpi_comm);  /* set to 0 x0 */
	
        local_map0 = (s2hat_flt8 *)calloc( nstokes*map_size,sizeof( s2hat_flt8));
	
        for( iter=0; (iter < niter) && (eps > epsilon*eps_ref); iter++) {

           /* pk -> Y pk */
       
           s2hat_alm2map( plms, cpixelization, cscan, nlmax, nmmax, nmvals, mvals, one, nstokes, first_ring, last_ring,
   	 	          map_size, local_map0, nlmax, local_alm1, nplm, plm, nprocs, myrank, mpi_comm);

           ypk_norm = mres_weighted_norm( cpixelization, cscan, nstokes, first_ring, last_ring, local_w8ring, map_size, local_map0, mpi_comm);

           /* ypk_norm = mres_norm( map_size*nstokes, local_map0, mpi_comm)*sqrt(cpixelization.parea[0]);    * this works for equal-area pixelizations only - e.g., HEALPIX */
	   
	   alpha_k = eps*eps/ypk_norm/ypk_norm;  /* if( myrank == root) { fprintf( stdout, "alpha = %.4e eps = %.4e ypk = %.4e [%.4e] \n", alpha_k, eps, ypk_norm, tmp); fflush(stdout); }*/

	   /* x_k+1 */

           combine_alm_sets( one, nstokes, nlmax, nmvals, mvals, nlmax, done, &local_alm[ishft], alpha_k, local_alm1, mpi_comm);
	   
	   /* r_k+1 */
	 
           s2hat_map2alm( plms, cpixelization, cscan, nlmax, nmmax, nmvals, mvals, one, nstokes, first_ring, last_ring,
  	                 local_w8ring, map_size, local_map0, nlmax, local_alm2, nplm, plm, nprocs, myrank, mpi_comm);

	   /* tmp=hres_norm( one, nstokes, nlmax, nmvals, mvals, nlmax, local_alm0, mpi_comm); if( myrank == root) { fprintf(stdout," |p_k| = %.4e \n", tmp); fflush( stdout); } */
	   
	   /* for( i=0; i< nalms; i++) { local_alm0[i].re -= alpha_k*local_alm2[i].re; local_alm0[i].im -= alpha_k*local_alm2[i].im; } */

           combine_alm_sets( one, nstokes, nlmax, nmvals, mvals, nlmax, done, local_alm0, -alpha_k, local_alm2, mpi_comm);

           if( test) {

              cls = (s2hat_flt8 *)calloc( (nlmax+1)*nspec, sizeof( s2hat_flt8));
              collect_cls( nmaps, imap, nstokes, nlmax, nmvals, mvals, nlmax, local_alm, 
                           nspec, cls, myrank, nprocs, root, mpi_comm);

	      if( myrank == root) {
 	         sprintf(fname, "cg_cls_iter_%d_map_%d.bin", iter, imap);
	         fout = fopen( fname, "w"); fwrite( cls, sizeof( s2hat_flt8), (nlmax+1)*nspec, fout); fclose( fout);
	      }
	      
	      free(cls);
	   }
	 
         /* compute the error the same on all procs */

	 epsold = eps;
         eps=hres_norm( one, nstokes, nlmax, nmvals, mvals, nlmax, local_alm0, mpi_comm);

         if( myrank == root) { fprintf( stdout, "PCG: finishing iteration : %d - map [%d], precision %.4e \n", iter, imap, eps/eps_ref); fflush( stdout); }
	 
	 if( eps < epsilon*eps_ref) break;
	 
         beta_k = eps*eps/epsold/epsold;
	 
         /* p_k+1 */
       
	 /* for( i=0; i< nalms; i++) { local_alm1[i].re = local_alm0[i].re+beta_k*local_alm1[i].re; local_alm1[i].im = local_alm0[i].im+beta_k*local_alm1[i].im; } */

         combine_alm_sets( one, nstokes, nlmax, nmvals, mvals, nlmax, beta_k, local_alm1, done, local_alm0, mpi_comm);  /* output -> local_alm1 */
	 
       }  /* over iterations */

       free( local_map0);
       free( local_alm0); free( local_alm1); free( local_alm2);
  
       iter_out[imap] = iter;
       eps_out[imap]  = eps/eps_ref;

     }   /* over maps */
     
  } /* if iterations are performed */
     
  return( 1);

}

s2hat_int4 get_random_alms( s2hat_int4 **streams, s2hat_int4 lmax, s2hat_flt8 *cl, s2hat_int4 nstokes, s2hat_int4 nmvals, s2hat_int4 *mvals, s2hat_dcomplex *alms)

{
  double *revect, *imvect;
  int i, is, ishft, l, lmax1, m, nlm;

  lmax1 = lmax + 1;
  nlm = lmax1*nmvals;

  i=0;
  for( is=0; is<nstokes; is++) {

    for( m=0; m<nmvals; m++) {
       
        revect = (double *)calloc( lmax1, sizeof( double));
        sprng_gaussran( streams[m], lmax1, revect);

        imvect = (double *)calloc( lmax1, sizeof( double));
        sprng_gaussran( streams[m], lmax1, imvect);

	ishft = is*nlm;
	
        if( mvals[m] == 0) {

          for( l = 2; l < lmax1; l++, i++) {
              alms[ is*nlm + l].re = sqrt(cl[is*lmax1+l])*revect[i++];    // real
              alms[ is*nlm + l].im = 0.0;                                 // imaginary
	  }

        } else {

          for( l = 2; l < lmax1; l++, i++) {
	      alms[ is*nlm + m*lmax1 + l].re = sqrt( 0.5*cl[is*lmax1+l])*revect[i];      // real
	      alms[ is*nlm + m*lmax1 + l].im = sqrt( 0.5*cl[is*lmax1+l])*imvect[i++];    // imaginary
	    }
        }

     }
  }

  free( revect);
  free( imvect);

}


void sprng_gaussran( int *stream, int nvect, double *rvect)

{
  double fac, rsq, v1, v2, tvect[2];
  int i;

  for( i=0; i< nvect/2; i++) {

    do
    {
       tvect[0] = sprng( stream);
       tvect[1] = sprng( stream);

       v1 = 2.0*tvect[0] - 1.0;
       v2 = 2.0*tvect[1] - 1.0;
       rsq = v1*v1 + v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);

    fac = sqrt(-2.0*log(rsq)/rsq);

    rvect[ 2*i]=v1*fac;
    rvect[ 2*i+1]=v2*fac;
  }

  if( nvect%2)
  {

    do
    {
       tvect[0] = sprng( stream);
       tvect[1] = sprng( stream);

       v1 = 2.0*tvect[0] - 1.0;
       v2 = 2.0*tvect[1] - 1.0;
       rsq = v1*v1 + v2*v2;
    }
    while( rsq>=1.0 || rsq==0.0);

    fac=sqrt( -2.0*log(rsq)/rsq);

    rvect[nvect-1]=v1*fac;
  }

}
