#ifndef ALLVAR_H
 #include "allvars.h"
#endif

void fft_3d(FFT_GRID *p);
void flip_data(FFT_GRID *p);
void factor_data(FFT_GRID *p);
void grid_transfer(int flag, FFT_GRID *f);
void copy_to_fft_array(int i, int j, int k, FFT_GRID *f, double *fft_data_local);
void copy_from_fft_array(int i, int j, int k, FFT_GRID * f, double *fft_data_local);

void fft_of_the_fields();

void init_all();
void allocate_grids();
void jing_simulation_info(int flag);

double linear_growth_factor(double a);
double dimensionless_growth_factor(double a);
double factor_x(double a);
double hubble_a(double a);
double integ_kernel(double a, void * param);
double integ_kernel_factor(double a);

void load_jing_pos();
void load_jing_vel();
void unit_pos();
void unit_vel();
void part_jing_pos(long int i, float * pos);
void part_jing_vel(long int i, float * pos);

void density_assignment();
void assign_part(float *pos, float mass);
void init_density_field();
void mass_assign(int i, int j, int k, double mass);
void prepare_the_density_field();

FFT_GRID *alloc_fft_grids();
FFT_GRID *pfft(FFT_GRID *head, int i, int j, int k);
void free_3d_grids();
void alloc_power_bin_arrays();
void free_srh_grids();
void free_all();

void out_put_data(POW_BIN *p_bin, int flag);

void power_calculator(FFT_GRID *f1, FFT_GRID *f2, POW_BIN *p_bin);
void init_power_bin(double k_min, double log_k_dis, POW_BIN *p_bin);
void add_power_bin(POW_BIN *a, POW_BIN *b);

void alloc_srh_grid();
SRH_GRID * p_srh(int i, int j, int k);
void init_srh_grids();

long int nearest_particle_id(float *input_pos);
void srh_range_nearest(int i_srh, int j_srh, int k_srh, int range, float *input_pos, double *best_r2, long int *best_id);
void srh_grid_local_nearest(int i, int j, int k, float *input_pos, double *local_best_r2, long int *local_best_id);
double distance_square(float *input_pos, float *part_pos, float *periodic_vector);
double find_min_surface_dis(int i, int j, int k, float *pos);

void velocity_assignment();
void init_velocity_field();
void velocity_copy(int i, int j, int k, float *vel);

void warn_and_end(char *s);
void state(char *s);

void make_window_function_data();
void make_spline_of_window_function();
double window_function(double k);
void correlation_delta_theta();
void correlation_delta_delta();
void correlation_theta_theta();
void create_theta_k();

void velocity_power();
void create_decomposed_velocity_fields();
void create_v_delta();
void create_v_stoch();
void create_v_b();