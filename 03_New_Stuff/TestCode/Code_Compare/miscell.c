void q_to_E(q, Euler)
double  q[4], Euler[3] ;
{
    double  A_q[3][3] ;
    void    q_to_A() ; 

    q_to_A(q, A_q) ;

    Euler[1] = acos(A_q[2][2])/PI*180.0 ;
    Euler[0] = atan2(A_q[2][0],-A_q[2][1])/PI*180.0 ;
      if(A_q[2][0] < 0.0) Euler[0]=360.0+Euler[0] ;
    Euler[2] = atan2(A_q[0][2],A_q[1][2])/PI*180.0 ;
      if(A_q[0][2] < 0.0)  Euler[2]=360.0+Euler[2] ;
}

void q_to_A(q, A)   
double q[4], A[3][3] ;
{ 
    
    A[0][0] = q[0]*q[0]-q[1]*q[1]-q[2]*q[2]+q[3]*q[3] ;
    A[0][1] = 2*(q[0]*q[1]+q[2]*q[3]) ;
    A[0][2] = 2*(q[0]*q[2]-q[1]*q[3]) ;
    A[1][0] = 2*(q[0]*q[1]-q[2]*q[3]) ;
    A[1][1] = -q[0]*q[0]+q[1]*q[1]-q[2]*q[2]+q[3]*q[3] ;
    A[1][2] = 2*(q[1]*q[2]+q[0]*q[3]) ;
    A[2][0] = 2*(q[0]*q[2]+q[1]*q[3]) ;
    A[2][1] = 2*(q[1]*q[2]-q[0]*q[3]) ;
    A[2][2] = -q[0]*q[0]-q[1]*q[1]+q[2]*q[2]+q[3]*q[3] ;
}

void antisym3(A, vec) 
double  vec[3], A[3][3] ;
{ 
    A[0][0] = A[1][1] = A[2][2] =     0.0 ;
    A[0][1] =  vec[2] ; A[1][0] = -vec[2] ;
    A[0][2] = -vec[1] ; A[2][0] =  vec[1] ;
    A[1][2] =  vec[0] ; A[2][1] = -vec[0] ; 
}

void antisym4(A, vec) 
double  vec[3], A[4][4] ;
{ 
    A[0][0] = A[1][1] = A[2][2] = A[3][3] = 0.0 ;
    A[0][1] =  vec[2] ; A[1][0] = -vec[2] ;
    A[0][2] = -vec[1] ; A[2][0] =  vec[1] ;
    A[1][2] =  vec[0] ; A[2][1] = -vec[0] ; 
    A[0][3] =  vec[0] ; A[3][0] = -vec[0] ;
    A[1][3] =  vec[1] ; A[3][1] = -vec[1] ;
    A[2][3] =  vec[2] ; A[3][2] = -vec[2] ;
}

void qcomp(a, b, c)   
double a[4], b[4],   /* input : a x b */  
       c[4] ;        /* output        */
{
   c[0] = a[3]*b[0]+a[2]*b[1]-a[1]*b[2]+a[0]*b[3] ;
   c[1] = -a[2]*b[0]+a[3]*b[1]+a[0]*b[2]+a[1]*b[3] ;
   c[2] = a[1]*b[0]-a[0]*b[1]+a[3]*b[2]+a[2]*b[3] ;
   c[3] = -a[0]*b[0]-a[1]*b[1]-a[2]*b[2]+a[3]*b[3] ;
}

void SYMINV3(a, ifail)
int     *ifail ;
double  a[3][3]  ;
{
      int     i, j, k, k1, r[3] ;
      double  big, test, kp1, km1, p[3], q[3] ;

      *ifail = 0 ; /*  construct truth table */
      for(i=0;i<3;i++)  r[i] = 1 ; /* begin inversion */
      for(i=0;i<3;i++)  {        /*  search for pivot */
	 big=0. ;
	 for( j = 0 ; j < 3 ; j++)  {
	      test = fabs(a[j][j]) ;
	      if( (test-big) > 0 ) {
		   if (r[j] < 0 ) { *ifail = 1  ;  return ; }
		   if (r[j] > 0 ) { big = test ;  k=j ; }
	      }
	 }
      /*  preparation for elimination step */
	 r[k]=0 ;
	 q[k]=1./a[k][k] ;
	 p[k]=1. ;
	 a[k][k]=0.0 ;
	 k1 = k+1 ;
	 kp1=k+1 ;
	 km1=k1-1 ;
	 if (km1 < 0)  { *ifail = 1 ; return ; }
	 if (km1 > 0)
	     for (j=0 ; j<km1 ; j++) {
		 p[j] = a[j][k] ;
		 q[j] = a[j][k] * q[k] ;
		 if (r[j] < 0)  { *ifail = 1 ; return ; }
		 if (r[j] > 0)  q[j] = -q[j] ;
		 a[j][k] = 0. ;
	     }
	 if ( (k1-3) > 0 ) { *ifail = 1 ; return ; }
	 if ( (k1-3) < 0 )
	     for (j=kp1 ; j < 3 ; j++) {
		 p[j] = a[k][j] ;
		 if (r[j] < 0) { *ifail = 1 ; return ; }
		 if (r[j] == 0)  p[j] = -p[j] ;
		 q[j] = -a[k][j]*q[k] ;
		 a[k][j] = 0.0 ;
	 }
      /*  elimination proper   */
	  for (j=0 ; j<3 ; j++) for (k=j ; k<3 ; k++)
		 a[j][k] = a[j][k] + p[j]*q[k] ;
      }

      /*  replace lower half   */
      for (j=0 ; j<3 ; j++) for (k=j ; k<3 ; k++) a[k][j] = a[j][k] ;
}

