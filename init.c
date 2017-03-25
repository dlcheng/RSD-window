#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "allvars.h"
#include "proto.h"
#include "define.h"

void init_all()
{
/* load postion and velocity files */  
  jing_simulation_info(FLAG);
  load_jing_pos();
  load_jing_vel();

/* allocate memory for data analysis */
  allocate_grids();
  alloc_power_bin_arrays();

/* init the srh grids */
  alloc_srh_grid();  
  init_srh_grids();
}    /* end init_all */

void allocate_grids()
{
  p_delta = alloc_fft_grids();
  p_vel_x = alloc_fft_grids();
  p_vel_y = alloc_fft_grids();   
  p_vel_z = alloc_fft_grids(); 
  p_theta = alloc_fft_grids(); 

/*
  p_vel_d_x = alloc_fft_grids(); 
  p_vel_d_y = alloc_fft_grids();
  p_vel_d_z = alloc_fft_grids();

  p_vel_s_x = alloc_fft_grids();
  p_vel_s_y = alloc_fft_grids();
  p_vel_s_z = alloc_fft_grids();

  p_vel_b_x = alloc_fft_grids();
  p_vel_b_y = alloc_fft_grids();
  p_vel_b_z = alloc_fft_grids();
*/

} /* end allocate_grids */ 
   
void jing_simulation_info(int flag)
{
  FILE *fp;

  char file_test[256];
  
  if(FILE_NUM == 1) 
    sprintf(file_test, "%s%s%04d%s%04d", INPUT_FOLDER, "/pos", SNAPSHOT, ".", CUR_STEP);
  else
    sprintf(file_test, "%s%s%04d%s%04d%s%02d", INPUT_FOLDER, "/pos", SNAPSHOT, ".", CUR_STEP,".",1);
  
  fp = fopen(file_test, "rb");
  rewind(fp);
  fseek(fp, sizeof(int), SEEK_SET);
  fread(&np, sizeof(np), 1, fp);                   /* get the total number of particles */
  fclose(fp);

 /*  Constants for Scale free Jing simulations
     P = a_j^(2/(ns+3))
     a_j at init_redshit is 1
     Pi at init_redshift is 1
     P = Pi + dp * Nstep 
     For normal CDM simulations
     it is effectively to set ns = -1
  */

  ng = NG;  
  Omega_m = OMEGA_M;
  Omega_v = OMEGA_V;
  boxsize = BOXSIZE;                                 /* in unit of Mpc/h */
  cur_step = CUR_STEP;
  total_step = TOTAL_STEP;
  init_z = INITIAL_Z;

  if (flag == 0)
    ns = NS;
  else
    ns = -1;

  dp = (pow(1.0 + init_z, 2.0/(ns+3.0)) - 1.0)/total_step;
  double j_factor = 1.0 + init_z;                    /* a_j / a = j_factor */

  double pa = 1.0 + dp * cur_step;
  double aj = pow(pa, (ns+3.0)/2.0);
  redshift = aj / j_factor;                          /* the redshift of the simulation */
  redshift = 1.0/redshift - 1.0;
  double normal_hubble = hubble_a(1.0/(1.0+redshift));
  vfactor = 2.0 / (ns+3) * pa * normal_hubble * boxsize;

  grid_dis = boxsize / (double) ng;  
  cell_volume = pow(grid_dis, 3);
  total_grid_num = ng * ng * ng;
  total_mass =  (3e10 * MPCTOM /(8.0 * Pi * G0) / SUNTOKG) * pow(boxsize, 3) * Omega_m;
  each_part_mass = total_mass / (double) np;
  
  printf("Simulation: %04d_%04d\n",SNAPSHOT, CUR_STEP);
  printf("Omega_m = %.2f, Omega_v = %.2f\n", Omega_m, Omega_v);
  printf("Redshift = %.2f\n", redshift);
  printf("Boxsize = %.1f [Mpc/h]\n",boxsize);
  printf("np= %ld\n",np);
  printf("Mt = %.3e [M_sun/h], Mp = %.3e [M_sun/h]\n", total_mass, each_part_mass);
  printf("Vfactor = %.3e\n\n", vfactor);

}  /* end load_jing_simulation_info */

