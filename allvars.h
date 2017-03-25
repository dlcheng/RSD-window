#ifndef ALLVAR_H
#define ALLVAR_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_fft_complex.h>

typedef struct fft_grid  FFT_GRID;
typedef struct pow_bin POW_BIN;
typedef struct link_node LINK_NODE;
typedef struct srh_grid SRH_GRID;

struct fft_grid                        /* the grid structure for FFT */
{
  float Re;                            /* the real part of the fft grid */
  float Im;                            /* the imaginary part of the fft grid */
};

struct pow_bin
{ 
  double k;                            /* middle k value of the bin */
  double p;	                           /* the mean power of the bin */
  double var_p;                        /* the variance of the power in the bin */
  double err_p;                        /* the error of the mean power in the bin */
  unsigned int n;                      /* number of grids in the bin */
};	

struct link_node
{
  long int id;                         /* the id of the particle from 0 to np-1 */
  LINK_NODE * next;                    /* the pointer to the next node */
};

struct srh_grid
{
  LINK_NODE * last;                    /* the pointer to the last node */
  LINK_NODE * first;                   /* the pointer to the first node */
};


/* Simulation information */
extern long int np;                    /* the total number of simulation particles */
extern double Omega_m;                 /* Omega_m at z=0 */
extern double Omega_m_z;               /* Omega_m at the redshift of the output */
extern double Omega_v;                 /* Omega_v at z=0 */
extern double Omega_v_z;               /* Omega_v at the redshift of the output */
extern double init_z;                  /* the initial redshift of the simulation */
extern double dp;                      /* P = Pi + dp * Nstep */
extern double redshift;                /* redshift of the output */
extern double boxsize;                 /* boxsize of the simulation */
extern double ns;                      /* the primordial power index of scale free simulation */
extern int cur_step;                   /* the current step of the simulation */
extern int total_step;                 /* the total steps of the simulation */

/* Particle information */
extern float * j_pos;                  /* the position info of Jing simulation */
extern float * j_vel;                  /* the velocity info of Jing simulation */
extern char  * reader_pos;             /* char arrays of pos */
extern char  * reader_vel;             /* char arrays of vel */

/* Velocity factor */
extern double vfactor;

/* THE GLOBAL SYBOLS USED IN THE CALCULATION */

extern FFT_GRID  *p_delta;             /* before FFT, this is the real space density constrast field */
extern FFT_GRID  *p_vel_x;             /* before FFT, this is the real space grid velocity in x direction */
extern FFT_GRID  *p_vel_y;             /* ....................................................y.......... */
extern FFT_GRID  *p_vel_z;             /* ....................................................z..........  */
extern FFT_GRID  *p_theta;             /* .........., this is the convergence field of the velocity */

extern FFT_GRID *p_vel_d_x;            /* the density correlated irrotational velocity part in x direction */
extern FFT_GRID *p_vel_d_y;            /* the density correlated irrotational velocity part in y direction */
extern FFT_GRID *p_vel_d_z;            /* the density correlated irrotational velocity part in z direction */
extern FFT_GRID *p_vel_s_x;            /* the density uncorrelated irrotational velocity part in x direction */
extern FFT_GRID *p_vel_s_y;            /* the density uncorrelated irrotational velocity part in y direction */
extern FFT_GRID *p_vel_s_z;            /* the density uncorrelated irrotational velocity part in z direction */ 
extern FFT_GRID *p_vel_b_x;            /* the rotational velocity part in x direction */
extern FFT_GRID *p_vel_b_y;            /* the rotational velocity part in y direction */
extern FFT_GRID *p_vel_b_z;            /* the rotational velocity part in z direction */

extern POW_BIN *p_bin_delta_delta;     /* power spectrum of the density contrast */
extern POW_BIN *p_bin_delta_theta;     /* power spectrum of the density and theta field */
extern POW_BIN *p_bin_theta_theta;     /* power spectrum of theta with theta */
extern POW_BIN *p_bin_vb_vb;           /* power spectrum of rotaional velocity part */
extern POW_BIN *p_bin_vd_vd;           /* power spectrum of the density correlated irrotational velocity part */
extern POW_BIN *p_bin_vs_vs;           /* power spectrum of the density uncoorelated irrotational velocity part */
extern POW_BIN *p_bin_window;          /* the bin value of window function */

extern long int ng;                    /* grid number in each dimension */
extern double grid_dis;
extern double cell_volume;
extern double total_mass;
extern double total_grid_num;
extern double each_part_mass;

/* spline of the window function */
extern gsl_interp_accel * spacc_window;
extern gsl_spline *sp_window;

/* related to the search grid */
extern SRH_GRID * p_srh_grid_array;    /* global pointer to the search grid */
extern long int srh_grid_num;          /* the total number of search grids */
extern double srh_grid_dis;            /* the size of srh_grid */
extern LINK_NODE * p_link_node_array;  /* the array of the link node with np node */

#endif
