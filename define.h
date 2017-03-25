/* 
 * unit used here:
 * 1 mass = 1 soloar mass /h
 * 1 length = 1 Mpc/h
 * 1 velocity = 100 km/s
 * such that H = h/Mpc velocity unit.
 */

#define  N_thread             128                         /* the number of threads */
#define  NG                   512                         /* the dimension of the grid, should be power of 2, like 128, 256, 512...*/
#define  N_srh                512                         /* the grid size of the search grid, could be any interger */
#define  BIN_NUMBER           20

/* the input file properties */
#define  OUTPUT_FOLDER        "./Result/"
#define  INPUT_FOLDER         "/home/dlcheng/data/2710"
#define  SNAPSHOT             2710
#define  CUR_STEP             5000
#define  FILE_NUM             1                           /* the number of input files */
#define  FLAG                 1                           /* 0 for scale free, 1 for LCDM */

/* Jing simulation parameters */
#define  OMEGA_M              1
#define  OMEGA_V              0
#define  NS                   0                           /* this is only used for the scale free simulations */
#define  INITIAL_Z            1029.15
#define  TOTAL_STEP           5000
#define  BOXSIZE              1200                        /* in unit of Mpc/h */

/* ************************************ */
#define  Pi                   3.1415926535897932384626433832795028842   
#define  SUNTOKG              1.9891e30
#define  MPCTOM               3.08567758e22
#define  G0                   6.67384e-11