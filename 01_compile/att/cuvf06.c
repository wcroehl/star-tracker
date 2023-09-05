int  cuvf(t_new, t_old, t_elapse, qq, q, P, w, b, 
          cnt, zsum, znum, ccdstep)
int    cnt ;
long   *znum ;
double ccdstep, t_new, *t_old, t_elapse, w[3], b[3], 
       qq[4], q[4], P[6][6], zsum[3] ; 
{
  int    i, j, k, ifail, E1, E2, ic ;
  double phi, w_mag, std_dev, temp, coef, 
         n_hat[3], q_propa[4], q_new[4], dq[4], dx[6], z[3],  
         M[4][4], P_propa[6][6], P_pro[6][6], R[3][3], B[3][3], 
	 H[3][6], HP[3][6], HPHT[3][3], PHT[6][3],   
         K[6][3], KH[6][6], I_KH[6][6], P_new[6][6], 
	 Q_k[6][6], Q_k2[6][6],
	 transit[6][6], A_pro[3][3], W_anti[3][3], ts[3][3], 
	 s_a, s_r, s_s, s_u, s_l, a_, b_, c_, ksai, eta, 
	 n33[3][3], n33sq[3][3], s_n33[3][3], c_n33sq[3][3],
	 del_t, ddel_t, Omega[3][3], OOmega[3][3],
	 Eq[4][3], Eqq[3], Edx[4], q_sum, W_pro[3],
	 w33[3][3], w_dg[3][3], w_up[3][3], w_lo[3][3], KR[6][6] ; 
  void   antisym3(), antisym4(), qcomp(), SYMINV3(), q_to_A() ;

  *znum = 0 ;
  for(i=0;i<3;i++) zsum[i] = 0.0 ;

  s_a = SIG_ARW*SIG_ARW ;
  s_r = SIG_RRW*SIG_RRW ;

  del_t = t_new - *t_old ;
  ddel_t = del_t*del_t ;
 
  for(i=0;i<3;i++) for(j=0;j<3;j++) Omega[i][j] = 0.0 ;
  for(i=0;i<3;i++) Omega[i][i] = w[i] ;

  for(i=0;i<3;i++) for(j=0;j<3;j++) OOmega[i][j] = 0.0 ;
  for(i=0;i<3;i++) for(j=0;j<3;j++) for(k=0;k<3;k++)
         OOmega[i][j] += Omega[i][k]*Omega[k][j] ; /* ??? */

  for(ic=0;ic<cnt;ic++) {
    /* for IST
    if(crf[ic]->mag >= 2.0) 
        std_dev = (4.5+0.7*(crf[ic]->mag-2.0)+ERR_POS)/3600.*PI/180.0 ;
    if(crf[ic]->mag <  2.0) 
    */

    std_dev = (6.5+ERR_POS)/3600.*PI/180.0 ;

    for(i=0;i<3;i++) for(j=0;j<3;j++) R[i][j] = 0.0 ;
    for(i=0;i<3;i++) R[i][i] = std_dev*std_dev ;

    if(fabs(t_new - *t_old) > 0.0001 ) 
      {       
      w_mag = sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]) ;     /* |ang. vel.| */
      for (i=0;i<3;i++) n_hat[i] = w[i]/w_mag ;    /* axis of rotation */ 
      phi = w_mag * (t_new - *t_old) ;   
      antisym4(M, n_hat) ;                   /* M : propagation matrix */

      for(i=0;i<4;i++) for(j=0;j<4;j++) M[i][j] *= sin(phi/2) ;
      for(i=0;i<4;i++) M[i][i] += cos(phi/2) ;

      for(i=0;i<4;i++) q_propa[i] = 0.0 ;
      for(i=0;i<4;i++) for(j=0;j<4;j++) q_propa[i] += M[i][j]*q[j] ;

      q_sum = sqrt(q_propa[0]*q_propa[0]+q_propa[1]*q_propa[1]
            +      q_propa[2]*q_propa[2]+q_propa[3]*q_propa[3]) ;
      for(i=0;i<4;i++) q_propa[i] /= q_sum ;

      antisym3(w33, w) ; /* UNNECESSARY */

      /* State transition matrix */
      antisym3(n33, n_hat) ;
      for(i=0;i<3;i++) for(j=0;j<3;j++) n33sq[i][j] = 0.0 ;
      for(i=0;i<3;i++) for(j=0;j<3;j++) for(k=0;k<3;k++)
           n33sq[i][j] += n33[i][k]*n33[k][j] ;
      for(i=0;i<3;i++) for(j=0;j<3;j++)
   	   s_n33[i][j] = n33[i][j]*sin(phi) ;
      for(i=0;i<3;i++) for(j=0;j<3;j++)
	   c_n33sq[i][j] = n33sq[i][j] * (1-cos(phi)) ;
      for(i=0;i<3;i++) for(j=0;j<3;j++)
          transit[i][j] = s_n33[i][j] + c_n33sq[i][j] ;
      for(i=0;i<3;i++) transit[i][i] += 1.0 ;

      a_ = (t_new - *t_old) ;   ksai = ( 1 - cos(phi) )/ phi ;
      b_ = a_ * ksai ;          eta  = 1 - sin(phi)/phi ;
      c_ = a_ * eta ;

      for(i=0;i<3;i++) for(j=3;j<6;j++)
      transit[i][j] = n33[i][j-3]*b_ + n33sq[i][j-3]*c_ ;
      for(i=0;i<3;i++) transit[i][i+3] += a_ ;

      for(i=0;i<3;i++) for(j=0;j<3;j++) ts[i][j] = transit[i][j+3] ; /* UNNECESSARY */

      for(i=3;i<6;i++) for(j=0;j<6;j++) transit[i][j] = 0.0 ;
      for(i=3;i<6;i++) transit[i][i] = 1.0 ;

      for(i=3;i<6;i++) for(j=0;j<6;j++) transit[i][j] = 0.0 ;
      for(i=3;i<6;i++) transit[i][i] = 1.0 ;

      /* P propagated */
      for(i=0;i<6;i++) for(j=0;j<6;j++) P_propa[i][j] = 0.0 ;
      for(i=0;i<6;i++) for(j=0;j<6;j++) for(k=0;k<6;k++)
  	      P_propa[i][j] += transit[i][k]*P[k][j] ;
      for(i=0;i<6;i++) for(j=0;j<6;j++)
             { P[i][j] = P_propa[i][j] ; P_propa[i][j] = 0.0 ; }
      for(i=0;i<6;i++) for(j=0;j<6;j++) for(k=0;k<6;k++)
              P_propa[i][j] += P[i][k] * transit[j][k] ;

      /* Q_k matrix */

      for(i=0;i<6;i++) for(j=0;j<6;j++) Q_k[i][j] = 0.0 ;

      /* 1 to 3 rows */
      for(i=0;i<3;i++)
                Q_k[i][i] = (s_a + s_r*ddel_t/3.0)*del_t ;
      for(i=0;i<3;i++)
                Q_k[i][i+3] = -s_r*ddel_t/2.0 ;

      /* 1 to 3 columns */
      for(i=0;i<3;i++)
                Q_k[i+3][i] = -s_r*ddel_t/2.0 ;

      /* diagonal terms */
      for(i=0;i<3;i++)
                Q_k[i+3][i+3] = s_r*del_t ;  /* ddel_t ? */

      for(i=0;i<6;i++) for(j=0;j<6;j++) P_propa[i][j] += Q_k[i][j] ; 
      
      for(i=0;i<6;i++) for(j=0;j<6;j++)
              P_pro[i][j] = (P_propa[i][j]+P_propa[j][i])/2.;

      /*      
      printf("\nP_pro\n") ;
      for(i=0;i<6;i++) {
         for(j=0;j<6;j++)
            printf(" %6.2e ",P_pro[i][j]) ;
	 printf("\n") ;
	 }
      printf("\n") ;
      */
      }  

  if(fabs(t_new - *t_old) < 0.0001 ) {
      for(i=0;i<6;i++) for(j=0;j<6;j++) P_pro[i][j] = P[i][j] ; 
      for(i=0;i<4;i++) q_propa[i] = q[i] ;
      }  
   
/* -----  update  ------ */

q_to_A(q_propa, A_pro) ;

for(i=0;i<3;i++) W_pro[i] = 0.0 ;
for(i=0;i<3;i++) for(j=0;j<3;j++) W_pro[i] += A_pro[i][j] * crf[ic]->L[j] ;
antisym3(W_anti, W_pro) ;

for(i=0;i<3;i++) for(j=0;j<3;j++) H[i][j] = W_anti[i][j] ;
for(i=0;i<3;i++) for(j=3;j<6;j++) H[i][j] = 0.0 ;

for(i=0;i<6;i++) for(j=0;j<3;j++) PHT[i][j] = 0.0 ;
for(i=0;i<6;i++) for(j=0;j<3;j++) for(k=0;k<6;k++)
        PHT[i][j] += P_pro[i][k]*H[j][k] ;

for(i=0;i<3;i++) for(j=0;j<3;j++)  HPHT[i][j] = 0.0 ;
for(i=0;i<3;i++) for(j=0;j<3;j++) for(k=0;k<6;k++)
        HPHT[i][j] += H[i][k]*PHT[k][j] ;

/* matrix B */
for(i=0;i<3;i++) for(j=0;j<3;j++) B[i][j] = HPHT[i][j] + R[i][j];
B[0][1] = B[1][0] = (B[0][1]+B[1][0])/2.0 ;
B[0][2] = B[2][0] = (B[0][2]+B[2][0])/2.0 ;
B[1][2] = B[2][1] = (B[1][2]+B[2][1])/2.0 ;

SYMINV3(B,&ifail) ;
if(ifail == 1) {
   printf(" Inverse of B matirx is singular\n") ;
   exit(1) ; }
B[0][1] = B[1][0] = (B[0][1]+B[1][0])/2.0 ;
B[0][2] = B[2][0] = (B[0][2]+B[2][0])/2.0 ;
B[1][2] = B[2][1] = (B[1][2]+B[2][1])/2.0 ;

for(i=0;i<6;i++) for(j=0;j<3;j++) K[i][j] = 0.0 ;
for(i=0;i<6;i++) for(j=0;j<3;j++) for(k=0;k<3;k++)
        K[i][j] += PHT[i][k]*B[k][j] ;

for(i=0;i<6;i++) for(j=0;j<6;j++) I_KH[i][j] = 0.0 ;
for(i=0;i<6;i++) I_KH[i][i] = 1.0 ;
for(i=0;i<6;i++) for(j=0;j<6;j++) for(k=0;k<3;k++)
		  I_KH[i][j] -= K[i][k]*H[k][j] ;

for(i=0;i<6;i++) for(j=0;j<6;j++) P[i][j] = 0.0 ;
for(i=0;i<6;i++) for(j=0;j<6;j++) for(k=0;k<6;k++)
        P[i][j] += I_KH[i][k]*P_pro[k][j] ;

for(i=0;i<6;i++) for(j=0;j<6;j++) P_new[i][j] = 0.0 ;
for(i=0;i<6;i++) for(j=0;j<6;j++) for(k=0;k<6;k++)
        P_new[i][j] += P[i][k]*I_KH[j][k] ;

for(i=0;i<6;i++) for(j=0;j<3;j++) KR[i][j] = 0.0 ;
for(i=0;i<6;i++) for(j=0;j<3;j++) for(k=0;k<3;k++)
        KR[i][j] += K[i][k]*R[k][j] ;
for(i=0;i<6;i++) for(j=0;j<6;j++) P_pro[i][j] = 0.0 ;
for(i=0;i<6;i++) for(j=0;j<6;j++) for(k=0;k<3;k++)
        P_pro[i][j] += KR[i][k]*K[j][k] ;
for(i=0;i<6;i++) if(P_pro[i][i] < 0.0 ) P_pro[i][i] = 0.0 ;
for(i=0;i<6;i++) for(j=0;j<6;j++)
        P_new[i][j] += P_pro[i][j] ;
for(i=0;i<6;i++) for(j=i;j<6;j++)
        P[i][j] = P[j][i] = (P_new[i][j] + P_new[j][i])/2. ;
for(i=0;i<6;i++) if(P[i][i] < 0.0 ) P[i][i] = 0.0 ;

/* dx update */
for(i=0;i<3;i++) z[i] = b_star[ic]->L[i] - W_pro[i] ;
for(i=0;i<3;i++) zsum[i] += sqrt(z[i]*z[i]) ;
                 (*znum)++ ;

for(i=0;i<6;i++) dx[i] = 0.0 ;
for(i=0;i<6;i++) for(j=0;j<3;j++) dx[i] += K[i][j]*z[j] ; 

temp = dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2] ; 
coef = 1. / sqrt(1.+ temp/4.) ;   
  
for(i=0;i<3;i++) dq[i] = -coef *dx[i]/2. ;  /* +coef */ 
                 dq[3] = coef ;
  
qcomp(dq,q_propa,q_new) ; 

q_sum = sqrt(q_new[0]*q_new[0]+q_new[1]*q_new[1]
      +      q_new[2]*q_new[2]+q_new[3]*q_new[3]) ;
for(i=0;i<4;i++) q_new[i] /= q_sum ;

for(i=0;i<3;i++) b[i] = b[i]+dx[i+3] ;

for(i=0;i<4;i++) q[i] = q_new[i] ;
*t_old = t_new ;

}

return(0) ;

} /* end of while loop */

