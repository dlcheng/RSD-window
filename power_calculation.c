#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "allvars.h"
#include "proto.h"
#include "define.h"

/* This routine calculates the self or cross correlation.
   the results are stored in p_bin.
*/
void power_calculator(FFT_GRID *f1, FFT_GRID *f2, POW_BIN *p_bin)
{
  int i, j, k;
  int i_n, j_n, k_n;
  double kx, ky, kz;
  int m;
  double loc_p;

#ifdef LOG_KBIN 
  double k_min = 2 * Pi / boxsize;
  double k_max = Pi / boxsize * ng * sqrt(3) / 2.0;		                    /* half the length of the three-dimension Nyquist frequency */
  double k_dis = (log10(k_max) - log10(k_min)) / (double) BIN_NUMBER;
#endif  

#ifdef LINEAR_KBIN
  double k_min = 0.;
  double k_max = 0.4;
  double k_dis = (k_max - k_min) / (double) BIN_NUMBER;
#endif  

  double k_now;
  double k_base = 2 * Pi / boxsize;
  
  init_power_bin(k_min, k_dis, p_bin);
  
  for(i = 0 ; i < ng; i++)
    for(j = 0; j < ng; j++)
      for(k = 0; k < ng; k++)
        {
		     i_n = i - ng/2;
		     j_n = j - ng/2;
		     k_n = k - ng/2;
		 	 
		     kx = k_base * i_n;
		     ky = k_base * j_n;
		     kz = k_base * k_n;	
		 		
	       k_now = sqrt(kx * kx + ky * ky + kz * kz);	     
           	     
	       if(k_now >= k_min && k_now <= k_max)
	         {
#ifdef LOG_KBIN            
			      m = (int) (log10(k_now/k_min) / k_dis);
#endif
#ifdef LINEAR_KBIN
            m = (int) ((k_now-k_min) / k_dis);
#endif            
			      if( m > BIN_NUMBER - 1)
			        m = BIN_NUMBER - 1;
			   
            loc_p = pfft(f1,i,j,k)->Re * pfft(f2,i,j,k)->Re + pfft(f1,i,j,k)->Im * pfft(f2,i,j,k)->Im;
		        p_bin[m].p += loc_p;
		        p_bin[m].n++;	   
           }			   
	    }
	    
  for(i = 0 ; i < BIN_NUMBER; i++)	    
	  p_bin[i].p = p_bin[i].p / (double) p_bin[i].n / pow(boxsize, 3);      /* divide the total volume */

/* do the survey of the k space another time to estimate the error */
  for(i = 0 ; i < ng; i++)
    for(j = 0; j < ng; j++)
      for(k = 0; k < ng; k++)
        {
         i_n = i - ng/2;
         j_n = j - ng/2;
         k_n = k - ng/2;
       
         kx = k_base * i_n;
         ky = k_base * j_n;
         kz = k_base * k_n; 
        
         k_now = sqrt(kx * kx + ky * ky + kz * kz);      
                 
         if(k_now >= k_min && k_now <= k_max)
           {
#ifdef LOG_KBIN            
            m = (int) (log10(k_now/k_min) / k_dis);
#endif
#ifdef LINEAR_KBIN
            m = (int) ((k_now-k_min) / k_dis);
#endif            
            if( m > BIN_NUMBER - 1)
              m = BIN_NUMBER - 1;
         
            loc_p = pfft(f1,i,j,k)->Re * pfft(f2,i,j,k)->Re + pfft(f1,i,j,k)->Im * pfft(f2,i,j,k)->Im;
            loc_p /= pow(boxsize, 3.);
            p_bin[m].var_p += pow(loc_p - p_bin[m].p, 2.);
           }         
      }

  for(i = 0 ; i < BIN_NUMBER; i++)
    {      
     p_bin[i].var_p = p_bin[i].var_p / ((double) p_bin[i].n - 1.);      /* the square variance */          
     p_bin[i].var_p = sqrt(p_bin[i].var_p);
     p_bin[i].err_p = p_bin[i].var_p / sqrt((double) p_bin[i].n);
    }
	
}    /* end collect_power_bin_data */	

void init_power_bin(double k_min, double k_dis, POW_BIN *p_bin)
{
  int i;
	
  for(i = 0 ; i < BIN_NUMBER; i++)
	{
#ifdef LOG_KBIN    
	 p_bin[i].k = pow(10., log10(k_min) + (i + 0.5) * k_dis);
#endif
#ifdef LINEAR_KBIN
   p_bin[i].k = k_min + (i + 0.5) * k_dis;  
#endif   
	 p_bin[i].p = 0;
	 p_bin[i].n = 0;
  }
}   /* end init_power_bin */	

/* add a to b and stored in b */
void add_power_bin(POW_BIN *a, POW_BIN *b)
{
  int i;

  for(i=0; i<BIN_NUMBER; i++)
    {
    b[i].p += a[i].p;
    b[i].n += a[i].n;  /* this information may not be useful */	
    }
}  /* end add_power_bin */