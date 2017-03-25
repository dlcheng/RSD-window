#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "allvars.h"
#include "proto.h"
#include "define.h"

/* the user interface is
   dvwf_j folder 2720(simulation) 5000(snapshot) ns 5000(total step number) 1200 (boxsize in unit of Mpc/h)
*/

int main(int argc, char *argv[])
{
  omp_set_num_threads(N_thread); 
  
  state("1. initialization");
  init_all();

  state("2. calculate density field");
  density_assignment();

  state("3. calculate velocity field");
  velocity_assignment();

  state("4. FFT of the fields");
  fft_of_the_fields();

  state("5. generate data of the window function");
  make_window_function_data();
  
  /*
  state("7. Decompose the velocity");
  velocity_power();
  */
  state("6. k, w(k), Delta^2(k)");
  out_put_data(p_bin_window, 1);
  /*
  state("9. Vd power k, P(k)");
  out_put_data(p_bin_vd_vd, 0);

  state("10. Vs power k, P(k)");
  out_put_data(p_bin_vs_vs, 0);

  state("11. Vb power k, P(k)");
  out_put_data(p_bin_vb_vb, 0);
  */
  free_all();

  return 1;
}   /* end main */
