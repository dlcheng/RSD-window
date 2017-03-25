#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "allvars.h"
#include "proto.h"
#include "define.h"


void density_assignment()
{
  long int i;
  float pos[3];

  init_density_field();

  for(i=0; i<np; i++)
   {
    part_jing_pos(i, pos);
    assign_part(pos, each_part_mass);
   }

  prepare_the_density_field();
  
}        /* end first_particle_survey */

void init_density_field()
{
  int i, j, k;

#pragma omp parallel shared(p_delta,ng) private(i,j,k)
{  
  #pragma omp for
  for(i=0; i<ng; i++)
    for(j=0; j<ng; j++)
      for(k=0; k<ng; k++)
      {
      pfft(p_delta,i,j,k)->Re = 0.0;
      pfft(p_delta,i,j,k)->Im = 0.0;
      }
}
} /* end init_density_field */  

void assign_part(float *pos, float mass)
{

  int i, j, k;

  i = (int) (pos[0] / grid_dis);
  if(i > ng-1)
    i = 0;

  j = (int) (pos[1] / grid_dis);
  if(j > ng-1)
    j = 0;

  k = (int) (pos[2] / grid_dis);
  if(k > ng-1)
    k = 0;

  mass_assign(i, j, k, mass);

}    /* end assign_part */

void mass_assign(int i, int j, int k, double mass)
{
  pfft(p_delta, i, j, k)->Re += mass;
}   /* end mass_assign */

void prepare_the_density_field()
{

  double each_grid_mean_mass = total_mass / total_grid_num;
  int i, j, k;
/* OPENMP */
#pragma omp parallel shared(ng, p_delta, each_grid_mean_mass) private(i, j, k)
{  
  #pragma omp for
  for(i=0; i<ng; i++)
    for(j=0; j<ng; j++)
       for(k=0; k<ng; k++)
         {
          pfft(p_delta,i,j,k)->Re = pfft(p_delta,i,j,k)->Re / each_grid_mean_mass - 1.0; /* how the density contrast is defined */
         }
}
} /* end prepare_the_density_field */