 #include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_fft_complex.h>

#include "allvars.h"

/* Simulation information */
 long int np;                    /* the total number of simulation particles */
 double Omega_m;                 /* Omega_m at z=0 */
 double Omega_m_z;               /* Omega_m at the redshift of the output */
 double Omega_v;                 /* Omega_v at z=0 */
 double Omega_v_z;               /* Omega_v at the redshift of the output */
 double init_z;                  /* the initial redshift of the simulation */
 double dp;                      /* P = Pi + dp * Nstep */
 double redshift;                /* redshift of the output */
 double boxsize;                 /* boxsize of the simulation */
 double ns;                      /* the primordial power index of scale free simulation */
 int cur_step;                   /* the current step of the simulation */
 int total_step;                 /* the total steps of the simulation */

/* Particle information */
 float * j_pos;                  /* the position info of Jing simulation */
 float * j_vel;                  /* the velocity info of Jing simulation */
 char *  reader_pos;             /* char arrays of pos */
 char *  reader_vel;             /* char arrays of vel */

/* Velocity factor */
 double vfactor;

/* THE GLOBAL SYBOLS USED IN THE CALCULATION */


 FFT_GRID  *p_delta;             /* before FFT, this is the real space density constrast field */
 FFT_GRID  *p_vel_x;             /* before FFT, this is the real space grid velocity in x direction */
 FFT_GRID  *p_vel_y;             /* ....................................................y.......... */
 FFT_GRID  *p_vel_z;             /* ....................................................z..........  */

 FFT_GRID  *p_theta;             /* .........., this is the convergence field of the velocity */

 FFT_GRID *p_vel_d_x;            /* the density correlated irrotational velocity part in x direction */
 FFT_GRID *p_vel_d_y;            /* the density correlated irrotational velocity part in y direction */
 FFT_GRID *p_vel_d_z;            /* the density correlated irrotational velocity part in z direction */
 FFT_GRID *p_vel_s_x;            /* the density uncorrelated irrotational velocity part in x direction */
 FFT_GRID *p_vel_s_y;            /* the density uncorrelated irrotational velocity part in y direction */
 FFT_GRID *p_vel_s_z;            /* the density uncorrelated irrotational velocity part in z direction */ 
 FFT_GRID *p_vel_b_x;            /* the rotational velocity part in x direction */
 FFT_GRID *p_vel_b_y;            /* the rotational velocity part in y direction */
 FFT_GRID *p_vel_b_z;            /* the rotational velocity part in z direction */


 POW_BIN *p_bin_delta_delta;     /* power spectrum of the density contrast */
 POW_BIN *p_bin_delta_theta;     /* power spectrum of the density and theta field */
 POW_BIN *p_bin_vb_vb;           /* power spectrum of rotaional velocity part */
 POW_BIN *p_bin_vd_vd;           /* power spectrum of the density correlated irrotational velocity part */
 POW_BIN *p_bin_vs_vs;           /* power spectrum of the density uncoorelated irrotational velocity part */
 POW_BIN *p_bin_window;          /* the bin value of window function */


 POW_BIN *p_bin_delta_delta;     /* power spectrum of the density contrast */
 POW_BIN *p_bin_delta_theta;     /* power spectrum of the density and theta field */
 POW_BIN *p_bin_theta_theta;     /* power spectrum of theta and theta */
 POW_BIN *p_bin_vb_vb;           /* power spectrum of rotaional velocity part */
 POW_BIN *p_bin_vd_vd;           /* power spectrum of the density correlated irrotational velocity part */
 POW_BIN *p_bin_vs_vs;           /* power spectrum of the density uncoorelated irrotational velocity part */
 POW_BIN *p_bin_window;          /* the bin value of window function */

 long int ng;                    /* grid number in each dimension */
 double grid_dis;
 double cell_volume;
 double total_mass;
 double total_grid_num;
 double each_part_mass;

/* spline of the window function */
 gsl_interp_accel * spacc_window;
 gsl_spline *sp_window;


/* related to the search grid */
 SRH_GRID * p_srh_grid_array;    /* global pointer to the search grid */
 long int srh_grid_num;          /* the total number of search grids */
 double srh_grid_dis;            /* the size of srh_grid */
 LINK_NODE *p_link_node_array;   /* the array of the link node with np node */
