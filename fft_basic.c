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
 * F is the real space real field, 
 * the FFT is to get the discrete fourier transform of the field, 
 * which is defined as 
 * F(k) = \int d^3 x F(x) exp(ik * x) over the periodic box of size L in each dimension.
 * The real k is mapped to the grid after 3d FFT as
 * k = (kg - ng/2) * 2\pi / L in each dimension, and the model value is
 * F(k) = V/N * F_{ng}(ng), where V is the volume of the box and N = ng^3 is the total number of grids.
 * Note that before apply the FFT, the F real space needs flipped.
*/

void fft_3d(FFT_GRID *p)
{  
  
  flip_data(p);
  grid_transfer(1, p);
  factor_data(p);
	
}     /* end fft_3d */	

/* the flip of the data is to make the grid (ng/2, ng/2, ng/2) effectively as k=0 of the Fourier mode */
void flip_data(FFT_GRID *p)
{
  int i, j, k;
  int m;
/* OPENMP */   
#pragma omp parallel shared(ng, p) private(i, j, k, m)
{ 
  #pragma omp for
  for(i = 0; i < ng; i++)
    for(j = 0; j < ng; j++)
      for(k = 0; k < ng; k++)
        {
		    m = i + j + k;	
         
	  	 if((m % 2) == 1)
		     {
		     pfft(p,i,j,k)->Re *= -1.0;
		     pfft(p,i,j,k)->Im *= -1.0;	
	       }	
	    }	
}      
}	        /* end flip_data */

/* the data is multiplied by the cell volume before estimating the fourier mode */
void factor_data(FFT_GRID *p)
{
  int i, j, k;  
/* OPENMP */
#pragma omp parallel shared(ng, p, cell_volume) private(i, j, k)
{  
  #pragma omp for
  for(i = 0; i<ng; i++)
    for(j = 0; j<ng; j++)
      for(k = 0; k < ng; k++)
         {
         pfft(p,i,j,k)->Re *= cell_volume;
         pfft(p,i,j,k)->Im *= cell_volume;
         }
}
}      /* end factor_data */

/* the worker who does the 3d FFT */
void grid_transfer(int flag, FFT_GRID *f)
{                             /* flag = 0 means k->x, 
                                 flag = 1 means x->k 
                                 f is the fft cube
                                 only flag = 1 is allowed here
                              */
 int i,j,k;
 double *fft_data_local;

#pragma omp parallel shared(ng, f, flag) private(i, j, k, fft_data_local)
{
 fft_data_local = (double *) malloc(2 * ng * sizeof(double));
 
 #pragma omp barrier
 #pragma omp for
 for(j=0; j<ng; j++)      
  {
   for(k=0; k<ng; k++)
     {
     copy_to_fft_array(-1,j,k,f,fft_data_local);

     if(flag == 0)
        gsl_fft_complex_radix2_inverse(fft_data_local,1,ng);
     if(flag == 1)	   
        gsl_fft_complex_radix2_forward(fft_data_local,1,ng);      

     copy_from_fft_array(-1,j,k,f,fft_data_local);     
     }
  }

  #pragma omp barrier 
  #pragma omp for
  for(k=0; k<ng; k++)
  {
   for(i=0; i<ng; i++)
     {
     copy_to_fft_array(i,-1,k,f,fft_data_local);

     if(flag == 0)
        gsl_fft_complex_radix2_inverse(fft_data_local,1,ng);
     if(flag == 1)
        gsl_fft_complex_radix2_forward(fft_data_local,1,ng);

     copy_from_fft_array(i,-1,k,f,fft_data_local);     
     }
  }

 #pragma omp barrier
 #pragma omp for
 for(i=0; i<ng; i++)
   {
   for(j=0; j<ng; j++)
     {
     copy_to_fft_array(i,j,-1,f,fft_data_local);

     if(flag == 0)
        gsl_fft_complex_radix2_inverse(fft_data_local,1,ng);
     if(flag == 1)
        gsl_fft_complex_radix2_forward(fft_data_local,1,ng);

     copy_from_fft_array(i,j,-1,f,fft_data_local);     
     }
    }
    
   free(fft_data_local);
}    
}                             /* end grid_transfer */

void copy_to_fft_array(int i, int j, int k, FFT_GRID *f, double *fft_data_local)
{
 int m;
 FFT_GRID *p_temp;

 if(i == -1)
   {
    for(m=0; m<ng; m++)
      {	  
      p_temp = pfft(f,m,j,k);
      fft_data_local[2*m] = p_temp->Re;
      fft_data_local[2*m+1] = p_temp->Im;
      }
   }

 if(j == -1)
   {
    for(m=0; m<ng; m++)
      {
      p_temp = pfft(f,i,m,k);
      fft_data_local[2*m] = p_temp->Re;
      fft_data_local[2*m+1] = p_temp->Im;
      }
   }

 if(k == -1)
   {
    for(m=0; m<ng; m++)
      {
      p_temp = pfft(f,i,j,m);
      fft_data_local[2*m] = p_temp->Re;
      fft_data_local[2*m+1] = p_temp->Im;
      }
   }
}                             /* end copy_to_fft_array */

void copy_from_fft_array(int i, int j, int k, FFT_GRID * f, double *fft_data_local)
{
 int m;
 FFT_GRID *p_temp; 
 if(i == -1)
   {
    for(m=0; m<ng; m++)
      {
      p_temp = pfft(f,m,j,k);        
      p_temp->Re = fft_data_local[2*m];
      p_temp->Im = fft_data_local[2*m+1];
      }
   }

 if(j == -1)
   {
    for(m=0; m<ng; m++)
      {
      p_temp = pfft(f,i,m,k);  
      p_temp->Re = fft_data_local[2*m];
      p_temp->Im = fft_data_local[2*m+1];
      }
   }

 if(k == -1)
   {
    for(m=0; m<ng; m++)
      {
      p_temp = pfft(f,i,j,m);  
      p_temp->Re = fft_data_local[2*m];
      p_temp->Im = fft_data_local[2*m+1];
      }
   }
}                             /* end copy_from_fft_array */