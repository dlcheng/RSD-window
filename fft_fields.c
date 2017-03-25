#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_fft_complex.h>

#include "allvars.h"
#include "proto.h"
#include "define.h"

void fft_of_the_fields()
{
/* Fourier tranform of the delta and velocity fields */	
  fft_3d(p_delta);
  fft_3d(p_vel_x);
  fft_3d(p_vel_y);
  fft_3d(p_vel_z);
}  /*end fft_of_the_fields */