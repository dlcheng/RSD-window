#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_fft_complex.h>

#include "allvars.h"
#include "proto.h"
#include "define.h"


/* the routines here generate the Fouier model of the decomposed velocity fields */

void velocity_power()
{
  create_decomposed_velocity_fields();
  POW_BIN *p_temp = (POW_BIN *) malloc(BIN_NUMBER * sizeof(POW_BIN));
  
  power_calculator(p_vel_d_x, p_vel_d_x, p_bin_vd_vd);
  power_calculator(p_vel_d_y, p_vel_d_y, p_temp);
  add_power_bin(p_temp, p_bin_vd_vd);
  power_calculator(p_vel_d_z, p_vel_d_z, p_temp);
  add_power_bin(p_temp, p_bin_vd_vd);    

  power_calculator(p_vel_s_x, p_vel_s_x, p_bin_vs_vs);
  power_calculator(p_vel_s_y, p_vel_s_y, p_temp);
  add_power_bin(p_temp, p_bin_vs_vs);
  power_calculator(p_vel_s_z, p_vel_s_z, p_temp);
  add_power_bin(p_temp, p_bin_vs_vs);    

  power_calculator(p_vel_b_x, p_vel_b_x, p_bin_vb_vb);
  power_calculator(p_vel_b_y, p_vel_b_y, p_temp);
  add_power_bin(p_temp, p_bin_vb_vb);
  power_calculator(p_vel_b_z, p_vel_b_z, p_temp);
  add_power_bin(p_temp, p_bin_vb_vb);    
  
  free(p_temp);
}  /* end velocity_power */

void create_decomposed_velocity_fields()
{
  create_v_delta();
  create_v_stoch();
  create_v_b();
}  /* end create_decomposed_velocity_fields */


void create_v_delta()
{
  int i, j, k;

  double kx, ky, kz;
  double k_base = 2 * Pi / boxsize;
  double k_now;
  double a = 1.0 / (1.0 + redshift);
  double normal_hubble_param = hubble_a(a);
  double f = dimensionless_growth_factor(a);
  double w_factor;

  double Real, Imag;

  for(i=0; i<ng; i++)
  	for(j=0; j<ng; j++)
  		for(k=0; k<ng; k++)
  		{
  		 kx = (i - ng/2) * k_base;
  		 ky = (j - ng/2) * k_base;
  		 kz = (k - ng/2) * k_base;	


         if (i==0 && j==0 && k==0)   /* skip k=0 mode */
         {
          pfft(p_vel_d_x,i,j,k)->Re = 0;
          pfft(p_vel_d_x,i,j,k)->Im = 0;

          pfft(p_vel_d_y,i,j,k)->Re = 0;
          pfft(p_vel_d_y,i,j,k)->Im = 0;

          pfft(p_vel_d_z,i,j,k)->Re = 0;
          pfft(p_vel_d_z,i,j,k)->Im = 0;

          continue;
         }

  		   k_now = sqrt(kx * kx + ky * ky + kz * kz);

         w_factor = window_function(k_now) * f;

         Real = pfft(p_delta,i,j,k)->Im * w_factor * normal_hubble_param / k_now / k_now * -1.0;
         Imag = pfft(p_delta,i,j,k)->Re * w_factor * normal_hubble_param / k_now / k_now;

         pfft(p_vel_d_x,i,j,k)->Re = Real * kx;
         pfft(p_vel_d_x,i,j,k)->Im = Imag * kx;

         pfft(p_vel_d_y,i,j,k)->Re = Real * ky;
         pfft(p_vel_d_y,i,j,k)->Im = Imag * ky;

         pfft(p_vel_d_z,i,j,k)->Re = Real * kz;
         pfft(p_vel_d_z,i,j,k)->Im = Imag * kz;

  		}

}  /* end create_v_delta */


void create_v_stoch()
{
  int i, j, k;
  double Real, Imag;
  double kx, ky, kz;
  double k_base = 2 * Pi / boxsize;
  double k_now;

  for(i=0; i<ng; i++)
  	for(j=0; j<ng; j++)
  		for(k=0; k<ng; k++)
  		{
         if(i==0 && j==0 && k==0)
         {
         pfft(p_vel_s_x,i,j,k)->Re = 0;
         pfft(p_vel_s_x,i,j,k)->Im = 0;

         pfft(p_vel_s_y,i,j,k)->Re = 0;
         pfft(p_vel_s_y,i,j,k)->Im = 0;

         pfft(p_vel_s_z,i,j,k)->Re = 0;
         pfft(p_vel_s_z,i,j,k)->Im = 0;

         continue;
         }
        
        kx = (i - ng/2) * k_base;
        ky = (j - ng/2) * k_base;
        kz = (k - ng/2) * k_base;
        
        k_now = sqrt(kx * kx + ky * ky + kz * kz);
        Real = pfft(p_vel_x,i,j,k)->Re * kx + pfft(p_vel_y,i,j,k)->Re * ky + pfft(p_vel_z,i,j,k)->Re * kz;
        Real = Real / k_now / k_now;

        Imag = pfft(p_vel_x,i,j,k)->Im * kx + pfft(p_vel_y,i,j,k)->Im * ky + pfft(p_vel_z,i,j,k)->Im * kz;
        Imag = Imag / k_now / k_now;

        pfft(p_vel_s_x,i,j,k)->Re = Real * kx - pfft(p_vel_d_x,i,j,k)->Re;
        pfft(p_vel_s_x,i,j,k)->Im = Imag * kx - pfft(p_vel_d_x,i,j,k)->Im;

        pfft(p_vel_s_y,i,j,k)->Re = Real * ky - pfft(p_vel_d_y,i,j,k)->Re;
        pfft(p_vel_s_y,i,j,k)->Im = Imag * ky - pfft(p_vel_d_y,i,j,k)->Im;

        pfft(p_vel_s_z,i,j,k)->Re = Real * kz - pfft(p_vel_d_z,i,j,k)->Re;
        pfft(p_vel_s_z,i,j,k)->Im = Imag * kz - pfft(p_vel_d_z,i,j,k)->Im;
  		}

}  		/* end create_v_s */

void create_v_b()
{
  int i, j, k;

  for(i=0; i<ng; i++)
  	for(j=0; j<ng; j++)
  		for(k=0; k<ng; k++)
  		{
  	    pfft(p_vel_b_x,i,j,k)->Re = pfft(p_vel_x,i,j,k)->Re - pfft(p_vel_d_x,i,j,k)->Re - pfft(p_vel_s_x,i,j,k)->Re;
  	    pfft(p_vel_b_x,i,j,k)->Im = pfft(p_vel_x,i,j,k)->Im - pfft(p_vel_d_x,i,j,k)->Im - pfft(p_vel_s_x,i,j,k)->Im;

  	    pfft(p_vel_b_y,i,j,k)->Re = pfft(p_vel_y,i,j,k)->Re - pfft(p_vel_d_y,i,j,k)->Re - pfft(p_vel_s_y,i,j,k)->Re;
  	    pfft(p_vel_b_y,i,j,k)->Im = pfft(p_vel_y,i,j,k)->Im - pfft(p_vel_d_y,i,j,k)->Im - pfft(p_vel_s_y,i,j,k)->Im;

  	    pfft(p_vel_b_z,i,j,k)->Re = pfft(p_vel_z,i,j,k)->Re - pfft(p_vel_d_z,i,j,k)->Re - pfft(p_vel_s_z,i,j,k)->Re;
  	    pfft(p_vel_b_z,i,j,k)->Re = pfft(p_vel_z,i,j,k)->Re - pfft(p_vel_d_z,i,j,k)->Re - pfft(p_vel_s_z,i,j,k)->Re;
  		}

}   /* end create_v_b */