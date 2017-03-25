#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "allvars.h"
#include "proto.h"
#include "define.h"
/* 
   the functions load the positions and velocities to the memory, only valid for big endian machines
   and return the particle structure float **j_pos and float **j_vel
*/

void load_jing_pos()
{
  FILE *fp;
  char file_name[256];
  long int readby;
  long int readby_total = 0L;
  int block_size_1, block_size_2, header_size;                                       /* this is 4 Bytes integer */
  long int file_len;
  int file_num= FILE_NUM;
  int i;
  
  if((reader_pos = (char *) malloc(np * 3 * sizeof(float))) == NULL)
    warn_and_end("Fail to allocate memory to pos.");

  for(i=1; i<=file_num; i++)
    {
    if(file_num == 1) 
      sprintf(file_name, "%s%s%04d%s%04d", INPUT_FOLDER, "/pos", SNAPSHOT, ".", CUR_STEP);
    else
      sprintf(file_name, "%s%s%04d%s%04d%s%02d", INPUT_FOLDER, "/pos", SNAPSHOT, ".", CUR_STEP,".",i);
   
    fp = fopen(file_name,  "rb");
    printf("Open %s\n", file_name);
    fseek(fp, 0, SEEK_END);
    file_len = ftell(fp);
    readby = 0L;   

    if(i==1)
      {
       header_size = 2*sizeof(long int) + 6 * sizeof(float) + 2*sizeof(int);
       fseek(fp, header_size, SEEK_SET);
       file_len -= header_size;
      }
    else
      rewind(fp);  
     
    while(readby < file_len)
       {
       fread(&block_size_1, sizeof(block_size_1), 1, fp);
       fread(&reader_pos[readby_total], sizeof(char), abs(block_size_1), fp);
       fread(&block_size_2, sizeof(block_size_2), 1, fp);
       
       if(abs(block_size_2) != abs(block_size_1))
         {
         printf("Inconsistency found in %s\n", file_name);
         warn_and_end("End");
         }

       readby += abs(block_size_1) + 2*sizeof(block_size_1);
       readby_total += abs(block_size_1);
       }  
    fclose(fp);  
    printf("Finish reading.\n");
    }  /* end for */

  j_pos = (float *) reader_pos;
  unit_pos();

  state("Done reading positions");
}   /* load the position file of Jing simulations */

void load_jing_vel()
{
  FILE *fp;
  char file_name[256];
  long int readby;
  long int readby_total = 0L;
  int block_size_1, block_size_2, header_size;                                       /* this is 4 Bytes integer */
  long int file_len;
  int file_num= FILE_NUM;
  int i;
  
  if((reader_vel = (char *) malloc(np * 3 * sizeof(float))) == NULL)
    warn_and_end("Fail to allocate memory to vel.");

  for(i=1; i<=file_num; i++)
    {
    if(file_num == 1) 
      sprintf(file_name, "%s%s%04d%s%04d", INPUT_FOLDER, "/vel", SNAPSHOT, ".", CUR_STEP);
    else
      sprintf(file_name, "%s%s%04d%s%04d%s%02d", INPUT_FOLDER, "/vel", SNAPSHOT, ".", CUR_STEP,".",i);
   
    fp = fopen(file_name,  "rb");
    printf("Open %s\n", file_name);
    fseek(fp, 0, SEEK_END);
    file_len = ftell(fp);
    readby = 0L;   

    if(i==1)
      {
       header_size = 2*sizeof(long int) + 6 * sizeof(float) + 2*sizeof(int);
       fseek(fp, header_size, SEEK_SET);
       file_len -= header_size;
      }
    else
      rewind(fp);  
     
    while(readby < file_len)
       {
       fread(&block_size_1, sizeof(block_size_1), 1, fp);
       fread(&reader_vel[readby_total], sizeof(char), abs(block_size_1), fp);
       fread(&block_size_2, sizeof(block_size_2), 1, fp);

       if(abs(block_size_2) != abs(block_size_1))
         {
          printf("Inconsistency found in %s\n", file_name);
          warn_and_end("End");
         }
         
       readby += abs(block_size_1) + 2*sizeof(block_size_1);
       readby_total += abs(block_size_1);
       }  
    fclose(fp);  
    printf("Finish reading.\n");
    }  /* end for */

  j_vel = (float *) reader_vel;
  unit_vel();

  state("Done reading velocities");
}  /* load the velocity file of Jing simulations */

void unit_pos()
{
 long int i;
#pragma omp parallel shared(np, j_pos, boxsize) private(i)
{
 #pragma omp for 
 for(i=0;i<3*np;i++)
   {
   j_pos[i] *= boxsize;
   }
}   
} /* end unit_pos */   

void unit_vel()
{
 long int i;
#pragma omp parallel shared(np, j_vel, vfactor) private(i)
{
 #pragma omp for 
 for(i=0;i<3*np;i++)
   {
   j_vel[i] *= vfactor;
   }
}   
} /* end unit_pos */   

void part_jing_pos(long int i, float * pos)
{
/* return the i-th particle's (starting from 0) position in pos */
  pos[0] = j_pos[i*3];
  pos[1] = j_pos[i*3 + 1];
  pos[2] = j_pos[i*3 + 2];
}       /* end part_jing_pos */ 

void part_jing_vel(long int i, float * vel)
{
/* return the i-th particle's (starting from 0) position in pos */
  vel[0] = j_vel[i*3];
  vel[1] = j_vel[i*3 + 1];
  vel[2] = j_vel[i*3 + 2];
}       /* end part_jing_pos */ 
