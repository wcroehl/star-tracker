int q_method(t, q, P, nstar, dist_index, out1a, out2) 
FILE   *out1a, *out2 ;
int    nstar ;
double t, q[4], P[3][3], *dist_index ;
{
    int    i, j, k, IFAIL ;      
    double V[Ns][3], W[Ns][3], B[3][3], S[3][3], A[3][3], 
           Z[3], AV[3], wgt[Ns], mag[Ns], 
           sig, lambda, std_dev, dev_tot ;
    void   q_to_A(), fastq(), SYMINV3() ;

    *dist_index = 0.0 ;
    for(i=0;i<nstar;i++) {
       for(j=0;j<3;j++) { V[i][j] = crf[i]->L[j] ; 
                          W[i][j] = b_star[i]->L[j] ; } ;
       mag[i] = crf[i]->mag ; }

    for(i=0;i<nstar;i++) wgt[i] = 1./nstar ;  

    /* distance index */
    for(i=0;i<nstar-1;i++) 
       *dist_index += 
         W[i][0]*W[i+1][0]+W[i][1]*W[i+1][1]+W[i][2]*W[i+1][2] ;
    *dist_index += 
         W[nstar-1][0]*W[0][0]+W[nstar-1][1]*W[0][1]+W[nstar-1][2]*W[0][2] ;

    for(i=0;i<3;i++)  for(j=0;j<3;j++) B[i][j] = 0.0 ;
    for(i=0;i<nstar;i++) for(j=0;j<3;j++) for(k=0;k<3;k++)  
             B[j][k] += wgt[i] * W[i][j] * V[i][k] ;

    sig = B[0][0] + B[1][1] + B[2][2] ; 
    for(i=0;i<3;i++) for(j=0;j<3;j++) S[i][j] = B[i][j] + B[j][i] ;

    Z[0] = Z[1] = Z[2] = 0.0 ;

    for(i=0;i<nstar;i++) {
       Z[0] += wgt[i]*(W[i][1]*V[i][2]-W[i][2]*V[i][1]);
       Z[1] += wgt[i]*(W[i][2]*V[i][0]-W[i][0]*V[i][2]);
       Z[2] += wgt[i]*(W[i][0]*V[i][1]-W[i][1]*V[i][0]); 
       }  
   
    fastq(S, Z, sig, &lambda, q) ;

    q_to_A(q, A) ;  

    for(i=0;i<3;i++) for(j=0;j<3;j++) P[i][j] = 0.0 ;
    for(i=0;i<nstar;i++) {
       for (j=0;j<3;j++) 
          AV[j] = A[j][0]*V[i][0]+A[j][1]*V[i][1]+A[j][2]*V[i][2] ;
       for (j=0;j<3;j++) for (k=0;k<3;k++)  
          P[j][k] += wgt[i]*AV[j]*AV[k];    
       } 

    for(i=0;i<3;i++) P[i][i] = 1 - P[i][i] ;  

    dev_tot = 0.0 ;
    for(i=0;i<nstar;i++) {
       if(mag[i] >= 2.0) std_dev = (4.5+0.7*(mag[i]-2.)+ERR_POS)/3600.*PI/180.;
       if(mag[i] <  2.0) std_dev = (4.5+ERR_POS)/3600.*PI/180. ; 
       dev_tot += 1./(std_dev*std_dev) ;
       }

    P[0][1] = P[1][0] = (P[0][1]+P[1][0])/2.0 ;
    P[0][2] = P[2][0] = (P[0][2]+P[2][0])/2.0 ;
    P[1][2] = P[2][1] = (P[1][2]+P[2][1])/2.0 ;
    SYMINV3(P,&IFAIL) ;
    if(IFAIL == 1)  { 
      fprintf(stderr, "Matrix inversion error at q_method()\n") ;
      return(10) ; }
    P[0][1] = P[1][0] = (P[0][1]+P[1][0])/2.0 ;
    P[0][2] = P[2][0] = (P[0][2]+P[2][0])/2.0 ;
    P[1][2] = P[2][1] = (P[1][2]+P[2][1])/2.0 ;

    for (i=0;i<3;i++) for (j=0;j<3;j++)
        P[i][j] /= dev_tot ; /* P_dthe_dthe = 4*P_dQ_dQ */
/*
    fprintf(out1a,"%15.6f",t) ;
    for(i=0;i<4;i++)
       fprintf(out1a,"  %15.10f",q[i]) ;
    fprintf(out1a,"\n") ;

    fprintf(out2,"%15.6f",t) ;
    fprintf(out2,"   %12.8f    %12.8f    %12.8f ", \
           sqrt(4*P[0][0])*180./PI*3600,  \
           sqrt(4*P[1][1])*180./PI*3600,  \
           sqrt(4*P[2][2])*180./PI*3600 ) ;    
    fprintf(out2,"\n") ;
*/
    return(0) ;
}

void fastq(S, Z, sig, lamb, X)  /* calculate optimal q */
double  S[3][3], Z[3], sig, *lamb, *X ; 
{
    int    i, j, k ;
    double kapa, delta, alpha, beta, xacc, q_denom,
           abcd[4], SZ[3], S2Z[3], S2[3][3], D[3][3],
           rtnewt() ;  
    void   funcd() ;

    kapa = S[0][0]*S[1][1]+S[0][0]*S[2][2]+S[1][1]*S[2][2] 
         - S[0][1]*S[1][0]-S[0][2]*S[2][0]-S[1][2]*S[2][1] ;

    delta = S[0][0]*(S[1][1]*S[2][2]-S[1][2]*S[2][1]) 
          + S[0][1]*(S[1][2]*S[2][0]-S[1][0]*S[2][2]) 
          + S[0][2]*(S[1][0]*S[2][1]-S[1][1]*S[2][0]) ;

    abcd[0] = sig*sig - kapa ;
    abcd[1] = sig*sig + Z[0]*Z[0]+Z[1]*Z[1]+Z[2]*Z[2] ;
    for(i=0;i<3;i++) SZ[i]=S[i][0]*Z[0]+S[i][1]*Z[1]+S[i][2]*Z[2] ;
    abcd[2] = delta + Z[0]*SZ[0]+Z[1]*SZ[1]+Z[2]*SZ[2] ;

    for(i=0;i<3;i++) S2Z[i]=S[i][0]*SZ[0]+S[i][1]*SZ[1]+S[i][2]*SZ[2] ;
    
    abcd[3] = Z[0]*S2Z[0]+Z[1]*S2Z[1]+Z[2]*S2Z[2] ;

    /* N-R Method to get the largest lambda with initial value  1  */
    xacc = 1.0e-8 ;
    *lamb = rtnewt(funcd,0.0,2.0,xacc,abcd,sig); /* X1,X2...interval */
    /* Calculate eigenvector corres. to lambda */
    alpha = *lamb**lamb - sig*sig + kapa ;
    beta = *lamb - sig ;
    X[3] = (*lamb+sig)*alpha - delta ;
    for (i=0;i<3;i++) for (j=0;j<3;j++) {
          D[i][j] = 0.0 ;
          if(i==j) D[i][j] = alpha ; }
   
    for(i=0;i<3;i++) for(j=0;j<3;j++) S2[i][j] = 0.0 ;
    for(i=0;i<3;i++) for(j=0;j<3;j++) for(k=0;k<3;k++)
          S2[i][j] += S[i][k]*S[k][j] ;
            
    for (i=0;i<3;i++) for (j=0;j<3;j++)
          S2[i][j] = beta * S[i][j] + D[i][j] + S2[i][j] ;

    for (i=0;i<3;i++) 
          X[i] = S2[i][0]*Z[0]+S2[i][1]*Z[1]+S2[i][2]*Z[2] ;

    q_denom = sqrt(X[0]*X[0]+X[1]*X[1]+X[2]*X[2]+X[3]*X[3]) ;

    for (i=0;i<4;i++) X[i] /= q_denom ;   
}

#define JMAX 20
double rtnewt(funcd, x1, x2, xacc, abcd, sig)
double x1, x2, xacc, abcd[4], sig ;
void (*funcd)() ;  
{ 
    int     j ;
    double  df, dx, f, rtn;

    rtn = 1.0 ; /* initial guess */
    for (j=1;j<=JMAX;j++) {
       (*funcd)(&f, &df, abcd, sig, rtn) ;
       dx   = f/df;
       rtn -= dx;
       if ((x1-rtn)*(rtn-x2) < 0.0)  
          fprintf(stderr, "Jumped out of brackets in RTNEWT\n") ; 
       if (fabs(dx) < xacc) return rtn ;  /* convergence */
    }
    fprintf(stderr, "Maximum number of iterations exceeded in RTNEWT\n") ;
}

void funcd(fn, df, abcd, sig, x)
double *fn, *df, abcd[4], sig, x ;
{
   double temp;
   
   temp = abcd[0]*abcd[1]+abcd[2]*sig-abcd[3] ;
   *fn  = x*x*x*x-(abcd[0]+abcd[1])*x*x-abcd[2]*x+temp ;
   *df  = 4*x*x*x-2*(abcd[0]+abcd[1])*x-abcd[2] ;
}


