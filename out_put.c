#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "allvars.h"
#include "proto.h"
#include "define.h"

void out_put_data(POW_BIN *p_bin, int flag)
{
/* flag = 0, only output to the terminal;
   flag = 1, also write a file */

  int i;
  char file_name[500];
  FILE *fp;

  if(flag == 1) 
   {
   sprintf(file_name, "%s%s%04d%s%04d%s%.2f%s%ld%s", \
           OUTPUT_FOLDER, \
           "wf",\
           SNAPSHOT, \
           "_", \
           CUR_STEP,\
           "_z_",\
           redshift, \
           "_", \
           ng, \
           ".txt");

   if((fp = fopen(file_name, "w+")) == NULL)    
     {
	   printf("Fail to open output file %s\n", file_name);
	   exit(0);
	   }
   }
  
  printf("%s\n","#(1)k   (2)P_dd   (3)P_dt   (4)P_tt" );  
  for(i = 0; i < BIN_NUMBER; i++)
    {			
      printf("%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n", \
                    p_bin[i].k, \
                    p_bin_delta_delta[i].p, \
                    p_bin_delta_theta[i].p, \
                    p_bin_theta_theta[i].p, \
                    p_bin_delta_delta[i].err_p, \
                    p_bin_delta_theta[i].err_p, \
                    p_bin_theta_theta[i].err_p);

      if(flag == 1)
        {    
        fprintf(fp, "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n", \
                    p_bin[i].k, \
                    p_bin_delta_delta[i].p, \
                    p_bin_delta_theta[i].p, \
                    p_bin_theta_theta[i].p, \
                    p_bin_delta_delta[i].err_p, \
                    p_bin_delta_theta[i].err_p, \
                    p_bin_theta_theta[i].err_p);
        }     
    }
  
  if(flag == 1)
   {
   fclose(fp);	
   }
}    /* end write file */	
