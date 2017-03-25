#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include "allvars.h"
#include "proto.h"
#include "define.h"

/* 
   the growth factor is got from O. Lahav et al. (1991) MNRAS 251, 128-136 
   f = dlnD/dlna works for any value of Omega_m and Omega_v.
*/

double linear_growth_factor(double a) /* the normalized linear growth factor */
{
  double x = factor_x(a);

  double integration_factor = integ_kernel_factor(a);

  double normal_factor = factor_x(1) * integ_kernel_factor(1);

  return (x / a * integration_factor) / normal_factor;

}  /* end linear_growth_factor */


double dimensionless_growth_factor(double a)
{
  double omega_m = Omega_m;
  double lambda_0 = Omega_v;
  double integration_factor = integ_kernel_factor(a);
  double x = factor_x(a);
  double y1 = lambda_0 * a * a - 0.5 * omega_m / a;

  double result;

  result = 1.0 / x / x * y1 - 1.0 + 1.0 / x / x / x * a / integration_factor;

  return result;
}  /* end dimensionless_growth_factor */


double hubble_a(double a)  /* this is the factor of H(a) / H(a=1) */
{
   double omega_m = Omega_m;
   double lambda_0 = Omega_v;

   double result;

   result = omega_m / pow(a, 3) + lambda_0 + (1 - omega_m - lambda_0) / a / a;

   return sqrt(result); 

}  /* end hubble_a */

double factor_x(double a)
{
  double omega_m = Omega_m;
  double lambda_0 = Omega_v;
  double result;

  result = (1.0 - omega_m - lambda_0) + omega_m / a + lambda_0 * a * a;
  return sqrt(result);                /* !!! notice the square root here !!!!*/
}  /* end factor_x */

double integ_kernel(double a, void * param)   
{
  double result = 1.0 / pow(factor_x(a) , 3);
  return result;
}    /* end integ_kenerl */

double integ_kernel_factor(double a)
{
  double result, error;
 
  gsl_integration_workspace * w = gsl_integration_workspace_alloc(10000);    
  gsl_function F;
  F.function = &integ_kernel;

  gsl_integration_qags(&F, 0, a, 0, 1e-4, 10000, w, &result, &error);     	   
  gsl_integration_workspace_free(w);

  return result;

}  /* end integ_factor */