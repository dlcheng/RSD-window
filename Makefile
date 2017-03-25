#--------------------------------------- Choice of K-bin
OPT += -DLINEAR_KBIN
#OPT += -DLOG_KBIN
#--------------------------------------- Select target computer
#SYSTYPE="Workstation"
SYSTYPE="Mac"
#SYSTYPE="uv2000"
#--------------------------------------- Adjust settings for target computer
ifeq ($(SYSTYPE),"Workstation")
CC       =   gcc   
OPTIMIZE =   -O3 -Wall 
GSL_INCL =  -I/home/dlcheng/Install/gsl/include
GSL_LIBS =  -L/home/dlcheng/Install/gsl/lib
endif

ifeq ($(SYSTYPE),"Mac")
CC       =   /usr/local/Cellar/gcc/5.2.0/bin/gcc-5  
OPTIMIZE =   -O3 -Wall
GSL_LIBS=   -L/usr/local/Cellar/gsl/1.16/lib 
GSL_INCL =  -I/usr/local/Cellar/gsl/1.16/include
endif

ifeq ($(SYSTYPE),"uv2000")
CC       =   gcc
OPTIMIZE =   -O3 -Wall
GSL_LIBS=   -L/home/dlcheng/lib/gsl-1.6/lib  -Wl
GSL_INCL =  -I/home/dlcheng/lib/gsl-1.6/include
endif

OPTIONS =  $(OPTIMIZE) $(OPT) 

EXEC   = dvwf_j
#        density-velocity-window-function for jing simulations

OBJS   = allvars.o fft_basic.o fft_fields.o init.o load_jing.o main.o linear_growth.o\
         mass_assign.o memory_control.o out_put.o power_calculation.o search_grid_init.o\
         search_grid_search.o velocity_assign.o velocity_decomp.o warn_end.o window_function.o

INCL   = allvars.h proto.h define.h Makefile


CFLAGS = $(OPTIONS) $(GSL_INCL) -fopenmp


LIBS   = $(GSL_LIBS) -lgsl -lgslcblas -lm -fopenmp

$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) 


clean:
	rm -f $(OBJS) $(EXEC) *.gch

