#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_fft_complex.h>

#include "allvars.h"
#include "proto.h"
#include "define.h"

/* 
 *  The rutines here calculate correlations of the fields, including
 *  1. the density power spectrum, stored in p_bin_delta_delta
 *  2. the density and theta power, stored in p_bin_delta_theta, 
 *  theta(k) = －ik * V(k) / H(a) / f， the Fourier transform has a different sign as the note,
 *  where H(a) has a redshift dependence.
 */

void make_window_function_data()
{
  correlation_delta_theta();
  correlation_delta_delta();
  correlation_theta_theta();

  int i;

  for(i=0; i<BIN_NUMBER; i++)
    {
    p_bin_window[i].k = p_bin_delta_delta[i].k;
    p_bin_window[i].p = p_bin_delta_theta[i].p / p_bin_delta_delta[i].p;  /* the window function is normalized by the dimensionless growth factor */
    p_bin_window[i].n = 0;                                                    /* this value is not used */
    }

  make_spline_of_window_function();                                           /* use spline function of window function */
}  /* end window_function_data */

void make_spline_of_window_function()
{
   double logk[BIN_NUMBER];
   double win[BIN_NUMBER];
   int i;

   for(i=0; i<BIN_NUMBER; i++)
    {
     logk[i] = log10(p_bin_window[i].k);
     win[i] = p_bin_window[i].p;	
    }

    spacc_window = gsl_interp_accel_alloc();
    sp_window = gsl_spline_alloc(gsl_interp_cspline, BIN_NUMBER);

    gsl_spline_init(sp_window, logk, win, BIN_NUMBER);

}  /* end make_spline_of_window_function */

double window_function(double k)
{                                                                           /* the User Interface of the window function, the input here is k*/
   double logk = log10(k);

   return gsl_spline_eval(sp_window, logk, spacc_window);
}  /* end window_function */

void correlation_delta_delta()
{
  power_calculator(p_delta, p_delta, p_bin_delta_delta);
}  /* end correlation_delta_delta */

void correlation_delta_theta()
{
  create_theta_k();
  power_calculator(p_delta, p_theta, p_bin_delta_theta);
} /* end correlation_delta_theta */

void correlation_theta_theta()
{
  power_calculator(p_theta, p_theta, p_bin_theta_theta);
}  /* end correlation_theta_theta */

void create_theta_k()
{
 int i, j, k;
 double kx, ky, kz;
 double k_base = 2.0 * Pi / boxsize;
 double a = 1.0 / (redshift + 1.0);
 double normal_hubble_param = hubble_a(a);
 double f = dimensionless_growth_factor(a);

 printf("f=%lf\n", f);
 
#pragma omp parallel shared(ng, k_base, p_theta, p_vel_x, p_vel_y, p_vel_z, normal_hubble_param, f) private(i,j,k,kx,ky,kz)
{
 #pragma omp for
 for(i=0; i<ng; i++)
 	for(j=0; j<ng; j++)
 	   for(k=0; k<ng; k++)
 	   {
       kx = k_base * (i - ng/2);
       ky = k_base * (j - ng/2);
       kz = k_base * (k - ng/2);
       
       pfft(p_theta,i,j,k)->Re = kx * pfft(p_vel_x,i,j,k)->Im + ky * pfft(p_vel_y,i,j,k)->Im + kz * pfft(p_vel_z,i,j,k)->Im;
       pfft(p_theta,i,j,k)->Im = -1.0 *(kx * pfft(p_vel_x,i,j,k)->Re + ky * pfft(p_vel_y,i,j,k)->Re + kz * pfft(p_vel_z,i,j,k)->Re);

       pfft(p_theta,i,j,k)->Re = pfft(p_theta,i,j,k)->Re / normal_hubble_param / f;         /* divide the time depdent hubble parameter */
       pfft(p_theta,i,j,k)->Im = pfft(p_theta,i,j,k)->Im / normal_hubble_param / f;
 	   }
}
}  /* end create_theta_k */