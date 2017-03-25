#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "allvars.h"
#include "proto.h"
#include "define.h"


void velocity_assignment()
{
  float input_pos[3];
  float vel[3];           /* the velocity of the nearest particle */
  long int n_id;          /* the nearest partilce ID */
  int i, j, k;  

  init_velocity_field();
/* OPENMP */
#pragma omp parallel shared(ng, grid_dis) private(i, j, k, input_pos, n_id, vel)
{  
  #pragma omp for  
  for(i=0; i<ng; i++)
  	for(j=0; j<ng; j++)
  	   for(k=0; k<ng; k++)
  	   	 {
          input_pos[0] = ((double) i + 0.5) * grid_dis;
          input_pos[1] = ((double) j + 0.5) * grid_dis;
          input_pos[2] = ((double) k + 0.5) * grid_dis;          
          n_id = nearest_particle_id(input_pos);
          part_jing_vel(n_id, vel);
          velocity_copy(i,j,k,vel);
  	   	 }
}  /* End OpenMP */
}  /* end velocity_assignment */

void init_velocity_field()
{
  int i, j,k;

#pragma omp parallel shared(ng, p_vel_x, p_vel_y, p_vel_z) private(i,j,k)
{
  #pragma omp for
  for(i=0;i<ng;i++)
    for(j=0;j<ng;j++)
      for(k=0;k<ng;k++)
      {
       pfft(p_vel_x,i,j,k)->Re = 0.0;
       pfft(p_vel_x,i,j,k)->Im = 0.0;

       pfft(p_vel_y,i,j,k)->Re = 0.0;
       pfft(p_vel_y,i,j,k)->Im = 0.0;

       pfft(p_vel_z,i,j,k)->Re = 0.0;
       pfft(p_vel_z,i,j,k)->Im = 0.0;  
      }
}
}  /* end init_velocity_field */

void velocity_copy(int i, int j, int k, float *vel)
{
  pfft(p_vel_x,i,j,k)->Re = vel[0];
  pfft(p_vel_y,i,j,k)->Re = vel[1];
  pfft(p_vel_z,i,j,k)->Re = vel[2];
} /* end velocity_copy */