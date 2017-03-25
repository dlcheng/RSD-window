#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "allvars.h"
#include "proto.h"
#include "define.h"


FFT_GRID * alloc_fft_grids()
{
  FFT_GRID *pointer;
  if((pointer = (FFT_GRID *)malloc(ng * ng * ng * sizeof(FFT_GRID))) == NULL)
     state("Fail to allocate memory of FFT grids."); 

  return pointer; 
}                             /* end alloc_fft_grids */

/* the the 1D array is mapped to the 3D array with index (i, j, k)
   note that the pointer to the grid is returned. */ 

FFT_GRID *pfft(FFT_GRID *head, int i, int j, int k)
{
  long int x_dis = ((long int) i) * ng * ng;
  long int y_dis = ((long int) j) * ng;
  long int z_dis = (long int ) k;
  long int total_dis = x_dis + y_dis + z_dis;

  return(&head[total_dis]);
}  /* end pfft */ 

void free_3d_grids()
{
  free(p_delta);
  free(p_vel_x);
  free(p_vel_y);
  free(p_vel_z);
  free(p_theta);

/*
  free(p_vel_d_x);
  free(p_vel_d_y);
  free(p_vel_d_z);

  free(p_vel_s_x);
  free(p_vel_s_y);
  free(p_vel_s_z);

  free(p_vel_b_x);
  free(p_vel_b_y);
  free(p_vel_b_z);
*/
}                               /* end free_3d_grids */
  
                            /* end alloc_fft_data */

void alloc_power_bin_arrays()
{
  p_bin_delta_delta = (POW_BIN *) malloc(BIN_NUMBER * sizeof(POW_BIN));
  p_bin_delta_theta = (POW_BIN *) malloc(BIN_NUMBER * sizeof(POW_BIN));
  p_bin_theta_theta = (POW_BIN *) malloc(BIN_NUMBER * sizeof(POW_BIN));
  p_bin_vb_vb = (POW_BIN *) malloc(BIN_NUMBER * sizeof(POW_BIN));
  p_bin_vd_vd = (POW_BIN *) malloc(BIN_NUMBER * sizeof(POW_BIN));
  p_bin_vs_vs = (POW_BIN *) malloc(BIN_NUMBER * sizeof(POW_BIN));
  p_bin_window = (POW_BIN *) malloc(BIN_NUMBER * sizeof(POW_BIN));
}	  /* end alloc_power_bin_arrays */

void free_srh_grids()
{ 
 free(p_srh_grid_array);
 free(p_link_node_array);
} /* end free_srh_grids */ 


void free_all()
{  
  free_3d_grids();

  free(reader_pos);
  free(reader_vel);
  
  free(p_bin_delta_delta);
  free(p_bin_delta_theta);
  free(p_bin_theta_theta);
  free(p_bin_vb_vb);
  free(p_bin_vd_vd);
  free(p_bin_vs_vs);
  free(p_bin_window);

  free_srh_grids();
		
}   /* end free_all */	
