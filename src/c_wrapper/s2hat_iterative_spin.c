
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"
#include "sys/types.h"
#include "sys/stat.h"
#include "sys/param.h"
#include "mpi.h"

#include "s2hat.h"


#define max( x, y) (((x) > (y)) ? (x) : (y))


s2hat_flt8 hres_norm( s2hat_int4 nmaps, s2hat_int4 nstokes, s2hat_int4 nlmax, s2hat_int4 nmvals, s2hat_int4 *mvals, s2hat_int4 lda, s2hat_dcomplex *local_alm, s2hat_flt8 *norm_result, MPI_Comm mpi_comm)

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
  *norm_result = rnorm;
  
  // return( rnorm);
  return 0  ;
}

s2hat_flt8 mres_norm( s2hat_int4 length, s2hat_flt8 *vect, s2hat_flt8 *norm_result, MPI_Comm mpi_comm)

/* computes a norm of a complex vector, cvect, of length elements *
 *                                               - rs, 2022/05/23 */
  
{
  s2hat_int4 i;
  s2hat_flt8 rnorm=0.0, rnorm_loc=0.0;

  for( i=0; i<length; i++) rnorm_loc += vect[i]*vect[i];

  MPI_Allreduce(&rnorm_loc, &rnorm, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

  rnorm = sqrt( rnorm);
  *norm_result = rnorm;

  // return( rnorm);
  return 0;
}

s2hat_flt8 mres_weighted_norm( s2hat_pixeltype cpixelization, s2hat_scandef cscan, s2hat_int4 nstokes, s2hat_int4 first_ring, s2hat_int4 last_ring, s2hat_flt8 *local_w8ring,
			       s2hat_int4 map_size, s2hat_flt8 *local_map, s2hat_flt8 *norm_result, MPI_Comm mpi_comm)

/* computes a weighted dot product of the distributed map         *
 * Tested only for HEALPIX                                        *
 *                                               - rs, 2022/06/02 */
  
{
  s2hat_int4 i, j, ipix, jpix, iring, is;
  s2hat_flt8 rnorm=0.0, rnorm_loc=0.0, rnorth, rsouth, tot_wght;

  // for( is = 0; is < nstokes; is++) {
  //   for( i = is*map_size, j = (is+1)*map_size-1, iring = first_ring; iring < last_ring+1; iring++) {
  //     tot_wght = cpixelization.parea[iring]*cpixelization.qwght[iring]*local_w8ring[iring-first_ring];

  //     if(cscan.nfl[iring]) {
  //        for( rnorth=0.0, ipix=0; ipix < cpixelization.nph[iring]; ipix++, i++) rnorth += local_map[i]*local_map[i];    /* north */
  //     } else
  //       i += cpixelization.nph[iring];
    
  //     if(cscan.sfl[iring]) {
  //        for( rsouth=0.0, jpix=0; jpix < cpixelization.nph[iring]; jpix++, j--) rsouth += local_map[j]*local_map[j];    /* south */
  //     } else
  //       j -= cpixelization.nph[iring];

  //     rnorm_loc += (rnorth+rsouth)*tot_wght;
  //   }
  // }
  for( is = 0; is < nstokes; is++) {
    for( i = is*map_size, j = (is+1)*map_size-1, iring = first_ring-1; iring < last_ring; iring++) {
      tot_wght = cpixelization.parea[iring]*cpixelization.qwght[iring]*local_w8ring[iring-first_ring+1];

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
  
  *norm_result = rnorm;
  // return( rnorm);
  return 0;
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

  return( 0);

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

      //  eps_ref = mres_norm( map_size*nstokes, &local_map[imap*nstokes*map_size], mpi_comm);
       mres_norm( map_size*nstokes, &local_map[imap*nstokes*map_size], &eps_ref, mpi_comm);
       
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

	  //  eps = mres_norm( map_size*nstokes, local_map0, mpi_comm);
     mres_norm( map_size*nstokes, local_map0, &eps, mpi_comm);
	   
           s2hat_map2alm( plms, cpixelization, cscan, nlmax, nmmax, nmvals, mvals, one, nstokes, first_ring, last_ring,
               	          local_w8ring, map_size, local_map0, nlmax, local_alm0, nplm, plm, nprocs, myrank, mpi_comm);

           /* update alms */

	   ishft = imap*nalms;
  	   /* for( i=0; i< nalms; i++) { local_alm[ishft+i].re -= local_alm0[i].re; local_alm[ishft+i].im -= local_alm0[i].im; }     this should be coded more precisely */

           combine_alm_sets( one, nstokes, nlmax, nmvals, mvals, nlmax, done, &local_alm[ishft], dmone, local_alm0, mpi_comm);
	   
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
     
  return( 0);

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

        // eps_ref=hres_norm( one, nstokes, nlmax, nmvals, mvals, nlmax, &local_alm[ishft], mpi_comm);
        hres_norm( one, nstokes, nlmax, nmvals, mvals, nlmax, &local_alm[ishft], &eps_ref, mpi_comm);

	// eps=hres_norm( one, nstokes, nlmax, nmvals, mvals, nlmax, local_alm0, mpi_comm);
  hres_norm( one, nstokes, nlmax, nmvals, mvals, nlmax, local_alm0, &eps, mpi_comm);

	// if( myrank == root) { fprintf( stdout, "PCG: initiating iteration : map [%d], precision %.4e \n", imap, eps/eps_ref); fflush( stdout); }
	
        local_alm1 = (s2hat_dcomplex *)calloc( nalms, sizeof( s2hat_dcomplex));    /* stores p_k */

        combine_alm_sets( one, nstokes, nlmax, nmvals, mvals, nlmax, dzero, local_alm1, done, local_alm0, mpi_comm);  /* copy r_0 to p_0 */
	
        local_alm2 = (s2hat_dcomplex *)calloc( nalms, sizeof( s2hat_dcomplex));   /* stores A p_k */
     
        for( iter=0; (iter < niter) && (eps > epsilon*eps_ref); iter++) {

           /* pk -> Y pk */
       
           s2hat_alm2map( plms, cpixelization, cscan, nlmax, nmmax, nmvals, mvals, one, nstokes, first_ring, last_ring,
   	 	          map_size, local_map0, nlmax, local_alm1, nplm, plm, nprocs, myrank, mpi_comm);

          //  ypk_norm = mres_weighted_norm( cpixelization, cscan, nstokes, first_ring, last_ring, local_w8ring, map_size, local_map0, mpi_comm);
           ypk_norm = mres_weighted_norm( cpixelization, cscan, nstokes, first_ring, last_ring, local_w8ring, map_size, local_map0, &ypk_norm, mpi_comm);

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

         /* compute the error the same on all procs */

	 epsold = eps;
        //  eps=hres_norm( one, nstokes, nlmax, nmvals, mvals, nlmax, local_alm0, mpi_comm);
         hres_norm( one, nstokes, nlmax, nmvals, mvals, nlmax, local_alm0, &eps, mpi_comm);

        //  if( myrank == root) { fprintf( stdout, "PCG: finishing iteration : %d - map [%d], precision %.4e \n", iter, imap, eps/eps_ref); fflush( stdout); }
	 
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
     
  return( 0);

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

        // eps_ref=hres_norm( one, nstokes, nlmax, nmvals, mvals, nlmax, &local_alm[ishft], mpi_comm);
        hres_norm( one, nstokes, nlmax, nmvals, mvals, nlmax, &local_alm[ishft], &eps_ref, mpi_comm);

	// eps=hres_norm( one, nstokes, nlmax, nmvals, mvals, nlmax, local_alm0, mpi_comm);
  hres_norm( one, nstokes, nlmax, nmvals, mvals, nlmax, local_alm0, &eps, mpi_comm);

	// if( myrank == root) { fprintf( stdout, "PCG: initiating iteration : map [%d], precision %.4e \n", imap, eps/eps_ref); fflush( stdout); }
	
    local_alm1 = (s2hat_dcomplex *)calloc( nalms, sizeof( s2hat_dcomplex));    /* stores p_k */
    combine_alm_sets( one, nstokes, nlmax, nmvals, mvals, nlmax, dzero, local_alm1, done, local_alm0, mpi_comm);  /* copy r_0 to p_0 */
    local_alm2 = (s2hat_dcomplex *)calloc( nalms, sizeof( s2hat_dcomplex));   /* stores A p_k */
    for( iter=0; (iter < niter) && (eps > epsilon*eps_ref); iter++) {
      /* pk -> Y pk */
      s2hat_alm2map_spin( cpixelization, cscan, spin, nlmax, nmmax, nmvals, mvals, one, first_ring, last_ring,
              map_size, local_map0, nlmax, local_alm1, nprocs, myrank, mpi_comm);

        //  ypk_norm = mres_weighted_norm( cpixelization, cscan, nstokes, first_ring, last_ring, local_w8ring, map_size, local_map0, mpi_comm);
          mres_weighted_norm( cpixelization, cscan, nstokes, first_ring, last_ring, local_w8ring, map_size, local_map0, &ypk_norm, mpi_comm);
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

        /* compute the error the same on all procs */

  epsold = eps;
      //  eps=hres_norm( one, nstokes, nlmax, nmvals, mvals, nlmax, local_alm0, mpi_comm);
        hres_norm( one, nstokes, nlmax, nmvals, mvals, nlmax, local_alm0, &eps, mpi_comm);
        // if( myrank == root) { fprintf( stdout, "PCG: finishing iteration : %d - map [%d], precision %.4e \n", iter, imap, eps/eps_ref); fflush( stdout); }

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
     
  return( 0);

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

        // eps = eps_ref = hres_norm( one, nstokes, nlmax, nmvals, mvals, nlmax, local_alm0, mpi_comm);
        hres_norm( one, nstokes, nlmax, nmvals, mvals, nlmax, local_alm0, &eps, mpi_comm);
        eps_ref = eps; 

	// if( myrank == root) { fprintf( stdout, "PCG: initiating iteration : map [%d], precision %.4e \n", imap, eps/eps_ref); fflush( stdout); }
	
        local_alm1 = (s2hat_dcomplex *)calloc( nalms, sizeof( s2hat_dcomplex));    /* stores p_k */

        combine_alm_sets( one, nstokes, nlmax, nmvals, mvals, nlmax, dzero, local_alm1, done, local_alm0, mpi_comm);  /* copy r_0 to p_0 */
	
        local_alm2 = (s2hat_dcomplex *)calloc( nalms, sizeof( s2hat_dcomplex));   /* stores A p_k */

        combine_alm_sets( one, nstokes, nlmax, nmvals, mvals, nlmax, dzero, &local_alm[ishft], done, local_alm2, mpi_comm);  /* set to 0 x0 */
	
        local_map0 = (s2hat_flt8 *)calloc( nstokes*map_size,sizeof( s2hat_flt8));
	
        for( iter=0; (iter < niter) && (eps > epsilon*eps_ref); iter++) {

           /* pk -> Y pk */
       
           s2hat_alm2map( plms, cpixelization, cscan, nlmax, nmmax, nmvals, mvals, one, nstokes, first_ring, last_ring,
   	 	          map_size, local_map0, nlmax, local_alm1, nplm, plm, nprocs, myrank, mpi_comm);

          //  ypk_norm = mres_weighted_norm( cpixelization, cscan, nstokes, first_ring, last_ring, local_w8ring, map_size, local_map0, mpi_comm);
           mres_weighted_norm( cpixelization, cscan, nstokes, first_ring, last_ring, local_w8ring, map_size, local_map0, &ypk_norm, mpi_comm);

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

	 
         /* compute the error the same on all procs */

	      epsold = eps;
        //  eps=hres_norm( one, nstokes, nlmax, nmvals, mvals, nlmax, local_alm0, mpi_comm);
         hres_norm( one, nstokes, nlmax, nmvals, mvals, nlmax, local_alm0, &eps, mpi_comm);

        //  if( myrank == root) { fprintf( stdout, "PCG: finishing iteration : %d - map [%d], precision %.4e \n", iter, imap, eps/eps_ref); fflush( stdout); }
	 
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
     
  return( 0);

}
