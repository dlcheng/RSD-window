#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "allvars.h"
#include "proto.h"
#include "define.h"

/*
 * this routine finds the nearest particle from the grid's chains 
 * the point is that if larger range is searched without finding nearer particle,
 * the routine is converged.
 * input_pos could be any postion  
 */

long int nearest_particle_id(float *input_pos) 
{
  int i_srh, j_srh, k_srh;  

  long int best_id;              
  double   best_r2; /* the theoretic largest distance is just 3/4 * L * L */
  long int next_range_best_id;
  double next_range_best_r2;
  int range = 0;

/* makes sure that input_pos is inbetween [0, boxsize) */
  while(input_pos[0]<0)
    input_pos[0] += boxsize;
  while(input_pos[0]>=boxsize)
    input_pos[0] -= boxsize;

  while(input_pos[1]<0)
    input_pos[1] += boxsize;
  while(input_pos[1]>=boxsize)
    input_pos[1] -=boxsize;

  while(input_pos[2]<0)
    input_pos[2] += boxsize;
  while(input_pos[2]>=boxsize)
    input_pos[2] -= boxsize;

/* which grid the input position belongs to */
  i_srh = (int) (input_pos[0]/srh_grid_dis);
  j_srh = (int) (input_pos[1]/srh_grid_dis);
  k_srh = (int) (input_pos[2]/srh_grid_dis);

 /* find the minimum distance to the suface of the grid */
  double min_dis = find_min_surface_dis(i_srh, j_srh, k_srh, input_pos); 

/* have a guess of the range 0 */
  srh_range_nearest(i_srh, j_srh, k_srh, range, input_pos,  &best_r2,  &best_id);

  if(best_r2 > min_dis * min_dis)        /* need to search the second lay */
    {
    while(1>0)
      {
      range++;
      srh_range_nearest(i_srh, j_srh, k_srh, range, input_pos, &next_range_best_r2, &next_range_best_id);
      
      if(next_range_best_r2 < best_r2)                 /* it is the best in the squared region */
         {
         best_r2 = next_range_best_r2;
         best_id = next_range_best_id;
         }  
      
      if(best_r2 <= pow(min_dis+range*srh_grid_dis,2)) /* it should also be best of the spherical region */
   	    break;
      }
     } 
  
  return best_id;  
}  /* end nearest_particle_id */

void srh_range_nearest(int i_srh, int j_srh, int k_srh, int range, float *input_pos, double *best_r2, long int *best_id)
{
  int i, j, k;
  double local_best_r2;
  long int local_best_id;
  *best_r2 = boxsize * boxsize;
  *best_id = -1;

  for(i = i_srh-range; i<=i_srh+range; i++)
  	for(j = j_srh-range; j<=j_srh+range; j++)
  	   for(k = k_srh-range; k<=k_srh+range; k++)
             {
             if(abs(i - i_srh) == range || abs(j - j_srh) == range || abs(k - k_srh) == range)
               {
               srh_grid_local_nearest(i, j, k, input_pos, &local_best_r2, &local_best_id);
               if(local_best_r2 < *best_r2)
              	{
                *best_r2 = local_best_r2;
                *best_id = local_best_id;		
               	}
               }
              }           

}  /* search the nearest particle at range, which starts from 0 */

void srh_grid_local_nearest(int i, int j, int k, float *input_pos, double *local_best_r2, long int *local_best_id)
{
/* the i, j, k here is not corrected for the periodic boundary condition 
   the grid here is guaranteed to contain at least one particle */

   int i_p, j_p, k_p;
   double local_r2;

   float periodic_vector[3];
   SRH_GRID * local_grid;
   LINK_NODE * local_node;
   float local_pos[3];

   i_p = i;
   j_p = j;
   k_p = k;

   while(i_p < 0)
   	 i_p += N_srh - 1;
   while(i_p > N_srh - 1)
     i_p -= N_srh -1;

   while(j_p < 0)
     j_p += N_srh - 1;
   while(j_p > N_srh - 1)
     j_p -= N_srh - 1;

   while(k_p < 0)
     k_p += N_srh - 1;
   while(k_p > N_srh -1)
     k_p -= N_srh -1;

   periodic_vector[0] = (i - i_p) * srh_grid_dis;
   periodic_vector[1] = (j - j_p) * srh_grid_dis;
   periodic_vector[2] = (k - k_p) * srh_grid_dis;

   local_grid = p_srh(i_p, j_p, k_p);
   local_node = local_grid->first;

/* incase the grid is empty */
   *local_best_r2 = boxsize * boxsize;
   *local_best_id = -1;

   while(local_node != NULL) /* this is a valid node */
     {
     part_jing_pos(local_node->id, local_pos); 	
     local_r2 = distance_square(input_pos, local_pos, periodic_vector);
     if(local_r2 < *local_best_r2)
       {
        *local_best_r2 = local_r2;
        *local_best_id = local_node->id;
       }
     local_node = local_node->next;
     }
}   /* this retruns the local nearest distance to the given poisiton */

double distance_square(float *input_pos, float *part_pos, float *periodic_vector)
{
  double r2=0;
  int i;

  for(i=0; i<3; i++)
  	r2 += pow(input_pos[i] - part_pos[i] - periodic_vector[i], 2);

  return r2;
}       /* end distance_square */

double find_min_surface_dis(int i, int j, int k, float *pos)
{
/* pos is inbetween [0, boxsize) */
  double min_dis = pos[0] - (double) i * srh_grid_dis;
  double dis_temp = ((double) i + 1.0) * srh_grid_dis - pos[0];

  if(dis_temp < min_dis)
     min_dis = dis_temp;

  dis_temp = pos[1] - (double) j * srh_grid_dis;
  if(dis_temp < min_dis)
     min_dis = dis_temp;

  dis_temp = ((double) j + 1.0) * srh_grid_dis - pos[1];
  if(dis_temp < min_dis)
     min_dis = dis_temp;

  dis_temp = pos[2] - (double) k * srh_grid_dis;
  if(dis_temp < min_dis)
     min_dis = dis_temp;

  dis_temp = ((double) k + 1.0) * srh_grid_dis - pos[2];
  if(dis_temp < min_dis)
     min_dis = dis_temp;  

  return min_dis;
} /* end find_min_surface_dis */ 