#include <time.h> 
#include <math.h> 
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <string.h>

#include "./nrutil.h"

#define  AU2KM  1.4959965e8 
#define  T   5670.9         
#define  Q_INI       1      /* 1: read from a file 2: input is given here */
#define  Ns   30            /* number of maximum stars measured in FOV */ 
#define  KF   -3600.0       
#define  Ctime     2.0      /* 5.0, time to change pm(), ckg() to dmt(),cuvf() */ 
#define  NUMSTAR   1854     /* 978 for mag 4.5, 1385 for mag 5.0 */
#define  Numg 2000
#define  TIMELIMIT 86390.0  

#define  SIG_ARW  1.0e-8   
#define  SIG_RRW  1.0e-8   

#define  TUNE1      10      
#define  TUNE2     1.0      

#define  TUNE_M4   150      
#define  TUNE_M2   300      

#define  DT_RMS    600  /* interval for the rms computation */
#define  T_STEP_CCD   0.1 
#define  T_STEP_GYRO  0.1 
#define  ERR_POS   0.3  /* median error for Hipa. is 0.77 mas-> set to 0.1 */
#define  ERR_MAG   0.5  /* 0.20+0.01, 0.35? */
#define  TOL6      50./3600.*3.141592653/180.0      /* 50 */
#define  TOL3      50./3600.*3.141592653/180.0      /* 50 */
#define  Limit     0.4  /* 1.0 magnitude limit : starID_pm() */
#define  DIST_TOL  360./3600.*3.141592653/180.0    /* for sym_test() */

#define  VICINITY  20      /* 5, 100   the vicinity for dmt() star ID  */
#define  mTOL0     0.8     /* magnitude limit0 for starID_dm */
#define  mTOL      0.8     /* magnitude limit  for starID_dm */
#define  BD_limit  6.0     /*        in starID_dm */
#define  FOV_limit 12.0    /*        in starID_dm */

#define  TOLERENCE 1.0e-9       /* iteration tolerence for batch.c */
#define  DT_LIMIT  200.0            /* batch interval */ 
#define  END_STAR  5671    /* 5792.0*15.  TIMELIMIT */

 double PI = 3.141592653589793100000 ;

 typedef struct {
       int  num_adj ;     /* the number of adjacent cells */
       int  *cell_num ;   /* adjacent cell number */
       }  ADJ_CELLS ;

 typedef struct {
       int  cell_id ;     /* the number of a cell */
       int  num_stars ;   /* the number of stars in the cell */
       int  *star_num ;   /* catalog number of stars */
       }  STAR_CELL ;

 typedef struct {
       int  cat_num ;     /* catalog number of a star */
       int  cat_mom ;   
       double ra ;       
       double dec ;
       double mag ;
       }  STAR ;

 struct list {
       int  starnum ;     /* the number of star for L[3] vector */
       double  IBL ;      /* the cosine of inter angle from dot prod */
       double  L[3] ;     /* this is only for generating test data */
       double  mag ;      /* magnitude(brightness) of the star */
       }  BLI[1000] ;

 struct VecStar {
       double  mag ;
       double  L[3] ;
       int     num ;
       int     id  ;
       } NEWBLI[Ns], NBLI, **crf, crf_temp[Ns] ;
              /* to store the star position vec. & mag. */

 struct VecStar2 {
       double  mag ;
       double  L[3] ;
       } **m_star, **b_star, **c_star ;
 
 STAR        stars[NUMSTAR] ;
 STAR_CELL   scell[225], scell2[34] ;
 ADJ_CELLS   adjcell[225]   ;

 /* Tracker(IST) to Body(GLAS)  
 double  T_B[3][3] = { 0.0,   -1.0,    0.0,  \
                       1.0,    0.0,    0.0,  \
                       0.0,    0.0,    1.0 } ;   
 */

 /* Tracker(IST) to GYRO */
 double  TR2B[3][3] = { 0.0,   -1.0,    0.0,  \
                       1.0,    0.0,    0.0,  \
                       0.0,    0.0,    1.0 } ;   
         
 /* Tracker(IST) to Body(IST) here */
 double  T_B[3][3] = { 1.0,    0.0,    0.0,  \
                       0.0,    1.0,    0.0,  \
                       0.0,    0.0,    1.0 } ;   
 
 double  B_EST_INI[] = {0.0, 0.0, 0.0}  ;      
 
 double  Q_INITIAL[] = {0.0, 0.0, 0.0, 1.0} ;  

