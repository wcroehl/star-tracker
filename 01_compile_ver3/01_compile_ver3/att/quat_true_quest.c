#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define FNLEN   256

main(int argc, char *argv[])
{
   FILE   *qtrue, *qesti, *delta, *bigD, *nstar ;   
   int    i, opp_sign ;     
   char   cat[18], ca[4][17], cqt[16], cstk[4], cq[4][17],
          fn_dir[FNLEN], fn_dir1[FNLEN], fn_dir2[FNLEN], fn_tmp[FNLEN] ;
   double t_a, t_q ;
   double pi, q_esti[4], q_true[4] ;
   double delqa[4], delthea[4], temp ; 
   void   qcomp() ;

   strcpy (fn_dir1, argv[1]);
   strcpy (fn_dir2, argv[2]);
   
   strcpy(fn_tmp, fn_dir1) ; strcat(fn_tmp, "/true_atti_lrs.dat") ;   
   qesti   = fopen(fn_tmp,"r") ;

   strcpy(fn_tmp, fn_dir2) ; strcat(fn_tmp, "/ckg_1a.dat") ; 
   qtrue = fopen(fn_tmp,"r") ;

   strcpy(fn_tmp, fn_dir2) ; strcat(fn_tmp, "/del_quest.dat") ;
   delta = fopen(fn_tmp,"w") ;

   strcpy(fn_tmp, fn_dir2) ; strcat(fn_tmp, "/big_qs_delta.dat") ;
   bigD = fopen(fn_tmp,"w") ;

   pi = 2.0*asin(1.0) ;

   fscanf(qtrue," %s ",cat) ;
   fscanf(qtrue," %s %s %s %s", ca[0],ca[1],ca[2],ca[3]) ;
   fscanf(qtrue," %s ",cat) ;

   while(fscanf(qesti,"%s", cqt) != EOF) 
     {
     fscanf(qesti," %s %s %s %s %s", cstk,cq[0],cq[1],cq[2],cq[3]) ;
     t_q = atof(cqt) ;    
     for(i=0;i<4;i++) q_true[i] = atof(cq[i]) ; 

     if( fabs(t_q - atof(cat)) < 0.001 ) 
         { 
         fscanf(qtrue,"%s %s %s %s", ca[0],ca[1],ca[2],ca[3]) ;
         t_a = atof(cat) ; 
         for(i=0;i<4;i++)  q_esti[i] = atof(ca[i]); 

         opp_sign = 0 ;
         for (i=0;i<4;i++) if (q_esti[i]*q_true[i] < 0) opp_sign++ ;
         if (opp_sign >= 3)  { /* sign = 1 ; */
                               for(i=0;i<4;i++) q_esti[i] = -q_esti[i] ;
                             }

         for(i=0;i<3;i++)  q_true[i] = -q_true[i] ; 

         qcomp(q_esti, q_true, delqa) ;

         temp = sqrt(delqa[0]*delqa[0]+delqa[1]*delqa[1]
                       +delqa[2]*delqa[2]+delqa[3]*delqa[3]) ;
         for(i=0;i<3;i++)  delthea[i] = 2 * delqa[i]/temp ;
         fprintf(delta, "%10.2f  %15.6f  %15.6f  %15.6f \n",
                t_q, delthea[0]*180./pi*3600, 
                     delthea[1]*180./pi*3600, 
                     delthea[2]*180./pi*3600) ;  
         if (fabs(delthea[0])*180./pi*3600. > 20.0  || fabs(delthea[1])*180./pi*3600. > 20.0)
	     fprintf(bigD, "%10.2f  %15.6f  %15.6f  %15.6f \n",
                     t_q, delthea[0]*180./pi*3600,
	    	     delthea[1]*180./pi*3600,
		     delthea[2]*180./pi*3600) ;

         fscanf(qtrue," %s ", cat) ; 
         }
     }

   fclose(qtrue) ;
   fclose(qesti) ;
   fclose(delta) ;
   fclose(bigD) ;
   return(0) ;
}
 
void qcomp(a, b, c) 
double a[4], b[4] ; 
double c[4] ;     
{
   c[0] = a[3]*b[0]+a[2]*b[1]-a[1]*b[2]+a[0]*b[3] ;
   c[1] = -a[2]*b[0]+a[3]*b[1]+a[0]*b[2]+a[1]*b[3] ;
   c[2] = a[1]*b[0]-a[0]*b[1]+a[3]*b[2]+a[2]*b[3] ;
   c[3] = -a[0]*b[0]-a[1]*b[1]-a[2]*b[2]+a[3]*b[3] ;
}

