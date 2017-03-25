#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "allvars.h"
#include "proto.h"
#include "define.h"


void alloc_srh_grid()
{
/* N_srh is the degree of the search grid in each dimension */	
 srh_grid_num = N_srh * N_srh * N_srh;
 p_srh_grid_array = (SRH_GRID *) malloc(srh_grid_num * sizeof(SRH_GRID));
 p_link_node_array = (LINK_NODE *)malloc(np * sizeof(LINK_NODE));

/* the search grid master node is initialized, OPENMP */

 long int i;
/* OPENMP */
#pragma omp parallel shared(srh_grid_num, p_srh_grid_array) private(i)
{
 #pragma omp for
 for(i=0; i<srh_grid_num; i++)
 	{
     p_srh_grid_array[i].last = NULL;        
     p_srh_grid_array[i].first = NULL;
 	}
}  
} /* end alloc_search_grid */


SRH_GRID * p_srh(int i, int j, int k)
{
/* this function return the pointer to a search grid at i, j, k */
  SRH_GRID *p;

/* this determines how the 1D array is mapped to the 3D array */

  long int x_dis = i * N_srh * N_srh;
  long int y_dis = j * N_srh;
  long int z_dis = k;
  long int total_dis = x_dis + y_dis + z_dis;

  p = &(p_srh_grid_array[total_dis]);

  return(p);

} /* end pointer_srh */

void init_srh_grids()
{
  srh_grid_dis = boxsize / (double) N_srh;

  long int id;
  float pos[3];         /* position of the partilce */
  int i, j, k;  
  LINK_NODE *local_node;
  SRH_GRID  *local_grid;

  for(id=0L; id<np; id++)
  	{
     part_jing_pos(id, pos);   /* load the position of the particle */   

     i = (int) (pos[0]/srh_grid_dis);
     if(i > N_srh-1)
     	i =0;
     j = (int) (pos[1]/srh_grid_dis);
     if(j > N_srh-1)
     	j = 0;
     k = (int) (pos[2]/srh_grid_dis);
     if(k > N_srh-1)
     	k = 0;
     
     local_grid = p_srh(i, j, k);
     local_node = &p_link_node_array[id];
     local_node->id =id;
     local_node->next = NULL;

     if(local_grid->first == NULL)  /* the first node of grid */
       {
        local_grid->first = local_node;
        local_grid->last = local_node;
       }
     else
       {
        local_grid->last->next = local_node;
        local_grid->last = local_node;
       }  
  	}
}  /* end init_srh_grids */