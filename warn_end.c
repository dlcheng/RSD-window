#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "allvars.h"
#include "proto.h"
#include "define.h"

void warn_and_end(char *s)
{
 printf("%s\n",s);
 exit(0);
}                           /* end warn_and_end */

void state(char *s)
{
 printf("%s\n",s);
}                           /* end state */
