starID_dm(t, q0, nstar, id_count, cctime, w, crf_count, id_star, id_body, 
          outi, outr, outo, outmag, outdist, ixy, x, y, xmag, t_start, t_end)
FILE   *outi, *outr, *outo,  *outmag, *outdist, *ixy ;
int    *id_count, nstar,  crf_count, *id_star, *id_body ; 
double t, q0[4], cctime, w[3], x[Ns], y[Ns], xmag[Ns], t_start, t_end ;
{  
   int     i, rtn, in_order(), find(), dmt_write() ; 
   double  *obs_ra, *obs_dec, *obs_mag, vici, Mp[3][3] ;
   void    dmt_atti() ; 

   obs_ra  = (double *) malloc(Ns * sizeof(double) ) ;
   if ( obs_ra  == NULL) {
      fprintf(stdout, " FAIL : malloc obs_ra (starID_dm) \n") ; 
      exit(5) ; }
   obs_dec = (double *) malloc(Ns * sizeof(double) ) ;
   if ( obs_dec == NULL) {
      fprintf(stdout, " FAIL : malloc obs_dec (starID_dm) \n") ;
      exit(5) ; }
   obs_mag = (double *) malloc(Ns * sizeof(double) ) ;
   if ( obs_mag == NULL) {
      fprintf(stdout, " FAIL : malloc obs_mag (starID_dm) \n") ;
      exit(5) ; }

   if (t > t_start && t < t_end) 
      {
      fprintf(stdout, "DM: t = %15.6f  q0: %12.8f %12.8f %12.8f %12.8f\n", 
			 t, q0[0], q0[1], q0[2], q0[3]) ;
      fprintf(stdout, "DM: t = %15.6f  w: %12.8f %12.8f %12.8f\n", 
			 t, w[0], w[1], w[2]) ;
      }

   dmt_atti(w, q0, Mp) ;
   rtn = 0 ;
   rtn = in_order(crf_count, nstar, Mp, obs_ra, obs_dec, obs_mag, x, y, xmag) ;
   if (rtn != 1)
      {
      fprintf(stdout,"starID_dm():in_order() didn't finish normally.\n") ;
      exit(6) ; 
      } 
  
   if (t > t_start && t < t_end) 
      {
      fprintf(stdout, "DM: t = %15.6f  cc_time = %15.6f\n", 
			 t, cctime) ;
      fprintf(stdout, "DM: crf_count = %4d  nstar = %4d  *id_count = %4d\n", 
		      crf_count, nstar, *id_count) ;
      }

   if (t < cctime+5.0)   vici = VICINITY*10.0 ;  /* cctime+1.0,  VI*5 */
   else                  vici = VICINITY+10.0 ;  /* VI+10.0 */ 
   
   if (t > t_start && t < t_end) 
      fprintf(stdout, "DM: t = %15.6f  cc_time = %15.6f  vici = %8.2f\n", 
			 t, cctime, vici) ;
   
   rtn = find(t, crf_count, nstar, obs_ra, obs_dec, obs_mag, 
        &(*id_count), id_star, id_body, vici, Mp, outdist, t_start, t_end) ;
   if (rtn != 1)
      {
      fprintf(stdout,"starID_dm():find() didn't finish normally.\n") ;
      exit(6) ; 
      } 

   rtn = dmt_write(t, &(*id_count), id_star, id_body, 
                outmag, outo, outi, outr, ixy, x, y, xmag) ;
   if (rtn != 1)
      {
      fprintf(stdout,"starID_dm():dmt_write() didn't finish normally.\n") ;
      exit(6) ; 
      } 

   free ( (void *) obs_mag) ;
   free ( (void *) obs_ra) ;
   free ( (void *) obs_dec) ;

   return(1) ;
} 

void dmt_atti(w, q0, Mp)  /* propagate attitude using the gyro meas. */ 
double w[3], q0[4], Mp[3][3] ;
 {
 int     i, j, k ;
 double  w_mag, n_hat[3], cwt, swt, qp[4], Ap[3][3] ;

 w_mag = sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]) ;
 for(i=0;i<3;i++)  n_hat[i] = w[i]/w_mag ;

 cwt = cos(w_mag*T_STEP_GYRO/2.0) ;  /* Originally multiplied by 10 */
 swt = sin(w_mag*T_STEP_GYRO/2.0) ;  /* Originally multiplied by 10 */

 qp[0] = cwt*q0[0] + swt*(n_hat[2]*q0[1] - n_hat[1]*q0[2] + n_hat[0]*q0[3]) ;
 qp[1] = cwt*q0[1] + swt*(-n_hat[2]*q0[0] + n_hat[0]*q0[2] + n_hat[1]*q0[3]) ;
 qp[2] = cwt*q0[2] + swt*(n_hat[1]*q0[0] - n_hat[0]*q0[1] + n_hat[2]*q0[3]) ;
 qp[3] = cwt*q0[3] + swt*(-n_hat[0]*q0[0] - n_hat[1]*q0[1] - n_hat[2]*q0[2])  ;

 q_to_A(qp, Ap) ;

 for(i=0;i<3;i++) for(j=0;j<3;j++) Mp[i][j] = 0.0 ;
 for(i=0;i<3;i++) for(j=0;j<3;j++) for(k=0;k<3;k++)
         Mp[i][j] += T_B[k][i] * Ap[k][j] ; /* Mp = Ap for IST */
}

int in_order(crf_count, obs_count, Mp, obs_ra, obs_dec, obs_mag, x, y, xmag)
int    crf_count, obs_count ;
double Mp[3][3], *obs_ra, *obs_dec, *obs_mag, x[Ns], y[Ns], xmag[Ns] ;
 {
 int    i, j, k ;
 double temporary, *crf_ra, *crf_dec, *crf_mag, tmp_x, tmp_y, tmp_mag ;
 struct VecStar  *temp  ;
 struct VecStar2 *temp2 ;

 temp = (struct VecStar *) malloc(sizeof(struct VecStar)) ; 
 if ( temp == NULL) {
   fprintf(stdout, " FAIL : malloc temp (starID_dm:in_order) \n") ; 
   exit(6) ; }
 temp2 = (struct VecStar2 *) malloc(sizeof(struct VecStar2)) ; 
 if ( temp2 == NULL) {
   fprintf(stdout, " FAIL : malloc temp (starID_dm:in_order) \n") ; 
   exit(6) ; }
 crf_ra  = (double *) malloc((crf_count) * sizeof(double) ) ;
 if ( crf_ra == NULL) {
   fprintf(stdout, " FAIL : malloc crf_ra (starID_dm:in_order) \n") ; 
   exit(6) ; }
 crf_dec = (double *) malloc((crf_count) * sizeof(double) ) ;
 if ( crf_dec == NULL) {
   fprintf(stdout, " FAIL : malloc crf_dec (starID_dm:in_order) \n") ; 
   exit(6) ; }
 crf_mag = (double *) malloc((crf_count) * sizeof(double) ) ;
 if ( crf_mag == NULL) {
   fprintf(stdout, " FAIL : malloc crf_mag (starID_dm:in_order) \n") ;
   exit(6) ; }

 for(i=0;i<obs_count;i++) { /* compute (RA,DEC) of observed stars */
   for(j=0;j<3;j++)
      c_star[i]->L[j] = 0.0 ;
   for(j=0;j<3;j++) for(k=0;k<3;k++)
      c_star[i]->L[j] += Mp[k][j] * m_star[i]->L[k] ;
   c_star[i]->mag = m_star[i]->mag ;

   obs_ra[i] = atan2(c_star[i]->L[1],c_star[i]->L[0])*180.0/PI ;
                  if ( obs_ra[i] < 0.0 ) obs_ra[i] += 360.0 ;
   obs_dec[i] = asin(c_star[i]->L[2])*180.0/PI ;
   obs_mag[i] = c_star[i]->mag ;
   } /* the end of for(i<obs_count) */

 for(i=0;i<obs_count;i++) { 
   for(j=i+1;j<obs_count;j++)  /* smaller first */
      if(obs_ra[j] < obs_ra[i])  
        {
        temporary = obs_ra[i] ;
        obs_ra[i] = obs_ra[j] ;
        obs_ra[j] = temporary ;

        temporary = obs_dec[i] ;
        obs_dec[i] = obs_dec[j] ;
        obs_dec[j] = temporary ;

        temporary = obs_mag[i] ;    
        obs_mag[i] = obs_mag[j] ;
        obs_mag[j] = temporary ;

        temp2->mag  = c_star[i]->mag ; 
        for(k=0;k<3;k++) temp2->L[k] = c_star[i]->L[k] ; 

        c_star[i]->mag  = c_star[j]->mag ; 
        for(k=0;k<3;k++) c_star[i]->L[k] = c_star[j]->L[k] ; 

        c_star[j]->mag  = temp2->mag ; 
        for(k=0;k<3;k++) c_star[j]->L[k] = temp2->L[k] ; 

        temp2->mag  = m_star[i]->mag ;  
        for(k=0;k<3;k++) temp2->L[k] = m_star[i]->L[k] ; 
	tmp_x = x[i] ;  tmp_y = y[i] ; tmp_mag = xmag[i] ;

        m_star[i]->mag  = m_star[j]->mag ; 
        for(k=0;k<3;k++) m_star[i]->L[k] = m_star[j]->L[k] ; 
	x[i] = x[j] ;   y[i] = y[j] ; xmag[i] = xmag[j] ;

        m_star[j]->mag  = temp2->mag ; 
        for(k=0;k<3;k++) m_star[j]->L[k] = temp2->L[k] ; 
	x[j] = tmp_x ;  y[j] = tmp_y ; xmag[j] = tmp_mag ;
        }
   } 

 for(i=0;i<crf_count;i++)
   {
   crf_ra[i]  = atan2(crf[i]->L[1],crf[i]->L[0])*180.0/PI ;
      if ( crf_ra[i] < 0.0 ) crf_ra[i] += 360.0 ;
   crf_dec[i] = asin(crf[i]->L[2])*180.0/PI ;
   crf_mag[i] = crf[i]->mag ;
   }

 for(i=0;i<crf_count;i++) 
   { 
   for(j=i+1;j<crf_count;j++)  /* smaller IBL first */
      if(crf_ra[j] < crf_ra[i])  
        {
        temporary = crf_ra[i] ;
        crf_ra[i] = crf_ra[j] ;
        crf_ra[j] = temporary ;

        temporary = crf_dec[i] ;
        crf_dec[i] = crf_dec[j] ;
        crf_dec[j] = temporary ;

        temporary = crf_mag[i] ; 
        crf_mag[i] = crf_mag[j] ;
        crf_mag[j] = temporary ;

        temp->mag  = crf[i]->mag ; 
        for(k=0;k<3;k++) temp->L[k] = crf[i]->L[k] ; 
        temp->num  = crf[i]->num ; 

        crf[i]->mag  = crf[j]->mag ; 
        for(k=0;k<3;k++) crf[i]->L[k] = crf[j]->L[k] ; 
        crf[i]->num  = crf[j]->num ; 

        crf[j]->mag  = temp->mag ; 
        for(k=0;k<3;k++) crf[j]->L[k] = temp->L[k] ; 
        crf[j]->num  = temp->num ; 
        }
   }

 free( (void *) crf_mag) ;  
 free( (void *) crf_dec) ;
 free( (void *) crf_ra) ;
 free( (void *) temp) ;  
 free( (void *) temp2) ;  
 return(1) ;
 }

int dmt_write(t, id_count, id_star, id_body, 
              outmag, outo, outi, outr, ixy, x, y, xmag)
FILE   *outmag, *outi, *outr, *outo, *ixy ;
int    *id_count, *id_star, *id_body  ;
double t, x[Ns], y[Ns], xmag[Ns] ;
 {
 int       i, j, k, rtn=0 ;
 double    alpha, delta, pos[3], ang_dist ;

 for(i=0;i<(*id_count);i++)
   {
   alpha = stars[*(id_star+i)-1].ra  / 180 * PI ;   /* RA_in_radian   */
   delta = stars[*(id_star+i)-1].dec / 180 * PI ;   /* Dec in radian */

   crf[i]->L[0] = pos[0] = cos(delta) * cos(alpha) ; /* usually on x axis */
   crf[i]->L[1] = pos[1] = cos(delta) * sin(alpha) ; /* usually on y axis */
   crf[i]->L[2] = pos[2] = sin(delta) ;  /* usually on z axis */
   crf[i]->mag  = stars[*(id_star+i)-1].mag ;
   crf[i]->num  = *(id_star+i) ;

   ang_dist = crf[i]->L[0]*pos[0]+crf[i]->L[1]*pos[1]+crf[i]->L[2]*pos[2] ;
   if (ang_dist >  1.0)  ang_dist =  1.0 ;
   if (ang_dist < -1.0)  ang_dist = -1.0 ;
   ang_dist = acos(ang_dist)*180.0/PI*3600.0; 

   for(j=0;j<3;j++) crf[i]->L[j] = pos[j] ; 

   for(j=0;j<3;j++) b_star[i]->L[j] = 0.0 ;
   for(j=0;j<3;j++) for(k=0;k<3;k++)
         b_star[i]->L[j] += T_B[j][k] * m_star[*(id_body+i)]->L[k] ;
   b_star[i]->mag = m_star[*(id_body+i)]->mag ;

   fprintf(outmag," %15.3f  2  %6d  %6.2f %6.2f %8.4f\n", 
            t, crf[i]->num, crf[i]->mag, b_star[i]->mag, ang_dist) ;

   fprintf(outi,"%15.6f  %20.15f %20.15f %20.15f %10.2f       %10d\n",
      t, crf[i]->L[0], crf[i]->L[1], crf[i]->L[2], crf[i]->mag, 
      crf[i]->num) ;
   
   fprintf(outo,"%15.6f  %20.15f %20.15f %20.15f %10.2f\n",
      t, m_star[*(id_body+i)]->L[0], m_star[*(id_body+i)]->L[1],
      m_star[*(id_body+i)]->L[2], m_star[*(id_body+i)]->mag) ;
   
   fprintf(ixy, "%15.6f  %10.2f %10.2f %10.2f %10.2f  %10d\n",  
      t, x[*(id_body+i)], y[*(id_body+i)], xmag[*(id_body+i)], m_star[*(id_body+i)]->mag, crf[i]->num) ;

   fprintf(outr,"%15.6f  %20.15f %20.15f %20.15f %10.2f\n",
      t, b_star[i]->L[0], b_star[i]->L[1],
      b_star[i]->L[2], b_star[i]->mag) ;
   }
 return(1) ;
 }

int find(t, crf_count, obs_count, obs_ra, obs_dec, obs_mag, 
             id_count, id_star, id_body, vici, Mp, outdist, t_start, t_end)
FILE    *outdist ;
int     crf_count, obs_count, *id_count, *id_star, *id_body ;
double  t, *obs_ra, *obs_dec, *obs_mag, vici, Mp[3][3], t_start, t_end ;
 {
 int     i, j, k, count, count0, sig, rtn=0,
         *find_num, *find_num2, *iptr, chosen_num, 
         new_obs() ;
 int     multi_obs() ;
 double  vec_dot, arc_dist, mag_diff ;

 count = sig = 0 ;

 if (t > t_start && t < t_end) 
   fprintf(stdout, "FIND: obs_count=%3d  *id_count=%3d\n", obs_count, *id_count) ; 
 for(i=0;i<obs_count;i++)
  {
  find_num   = (int *) malloc(10*sizeof(int)) ; 
  if (find_num == NULL) { 
     fprintf(stdout, " FAIL : malloc find_num (starID_dm:find) \n") ; 
     exit(5); 
     }

 if (t > t_start && t < t_end) 
   fprintf(stdout, "Here1: *id_count=%3d\n", *id_count) ; 

  iptr = find_num ;

 if (t > t_start && t < t_end) 
  fprintf(stdout, "FIND: t=%15.6f c_star[%2d]->mag = %7.2f\n", t, i, c_star[i]->mag) ; 
  for(j=0;j<crf_count;j++)
    {
    for(k=0;k<(*id_count);k++) 
       if (crf[j]->num == id_star[k]) sig = 8 ; 
    if (sig == 8)  { sig = 0 ; continue ; } 

    vec_dot = c_star[i]->L[0]*crf[j]->L[0]
            + c_star[i]->L[1]*crf[j]->L[1]
            + c_star[i]->L[2]*crf[j]->L[2] ;

    if (vec_dot > 1.0)  vec_dot =  1.0 ;
    if (vec_dot < -1.0) vec_dot = -1.0 ;
    arc_dist = acos(vec_dot)*180.0/PI*3600.0 ;
    mag_diff = (c_star[i]->mag-crf[j]->mag) ;

    if (t > t_start && t < t_end) 
      {
      fprintf(stdout, "FIND: t=%15.6f  crf[%2d]->num=%5d  mag=%5.2f ", t, j, crf[j]->num, crf[i]->mag) ;
      fprintf(stdout, "arc_dist (arcsec) =%10.4f\n", arc_dist) ;
      }

    if(fabs(arc_dist) < vici && fabs(mag_diff) < mTOL0) /* found one candidate: was VICINITY */
      {
      *iptr = stars[crf[j]->num-1].cat_num ;
      fprintf(outdist, "1  %15.6f %5d %15.8f \n", t, *iptr, arc_dist) ;
      count++  ;
      iptr++ ;
      } 
    }  /* end of for(j) */

    switch(count)
      {
      case 0 : /* False observation or New observation */
               find_num2   = (int *) malloc(50*sizeof(int)) ;  
               if (find_num2 == NULL) { 
                  fprintf(stdout," FAIL : malloc find_num2 (starID_dm:find)\n");
                  exit(5); 
                  }
 
               count0 = 0 ;
               chosen_num = crf[(int) floor(crf_count/2.0)]->num ;

               rtn = new_obs(t, t_start, t_end, &count0, find_num2, chosen_num,
                    *(c_star+i), (obs_ra+i), (obs_dec+i), (obs_mag+i),  
                    id_star, *id_count, vici, Mp, outdist) ;

               if (rtn != 1)
                  {
                  fprintf(stdout,"starID_dm()::new_obs() finished abnormally.\n") ;
                  exit(6) ; 
                  } 

               switch(count0)
                 {
                 case 0 : sig = 0 ; /* No identified star is reported */
                          break ;
                 case 1 : sig = 10; /* New star : *find_num2  */
                          chosen_num = *find_num2 ;
                          break ;
                 default: sig = 20; /* Multiple candidates : to choose one */
                          chosen_num = multi_obs(t, *(c_star+i), find_num2, count0) ;  
                 } ;
               free ( (void *) find_num2) ;
               break ;

      case 1 : sig = 1 ; /* New star : *find_num */
               chosen_num  = *find_num ;
               break ;

      default: sig = 2 ; /* Multiple candidates : to choose one */
               chosen_num = multi_obs(t, *(c_star+i), find_num, count) ; 
      }

  free( (void *) find_num ) ; 

  if(sig != 0)  {  
                id_star[*id_count] = chosen_num ;
                id_body[(*id_count)++] = i ;  /* for c_star */
                }
  count = 0 ;
  } /* end of for(i) */

 return(1) ;
 }
   
int new_obs(t, t_start, t_end, count, candi_num, oldnum,
         new, new_ra, new_dec, new_mag, id_star, id_count, vici, Mp, outdist)
FILE   *outdist ;
struct VecStar *new ;
int    *count, *candi_num, oldnum, *id_star, id_count ;
double t, t_start, t_end, *new_ra, *new_dec, *new_mag, vici, Mp[3][3] ; 
 {
 static int  previ_N_zone, thisnum ;
 int    i, j, *intptr, N_zone, astaragain, thisnum0, thisnum1, number,
        zone() ;
 double BD_vec[3], BD_ra, BD_dec, 
        ra_diff_new, ra_diff_BD, i_ra, i_dec, i_mag, i_L[3],
        vec_dot, arc_dist ;

 astaragain = 0 ;
 for(i=0;i<3;i++) BD_vec[i] = Mp[2][i] ; /* give star vector in the CRF frame */ 
 BD_ra = atan2(BD_vec[1], BD_vec[0]) * 180.0 / PI ;
               if ( BD_ra  <  0.0) BD_ra += 360.0 ;
 BD_dec = asin(BD_vec[2]) * 180.0 / PI ;

 N_zone = zone(BD_dec) ; /* N_zone_ON */

 if (N_zone == -1)
   {
   fprintf(stdout,"find()::zone() didn't finish normally.\n") ;
   exit(6) ; 
   } 

 if (N_zone != previ_N_zone)
   {
   for(j=0 ; j < scell2[N_zone-1].num_stars ; j++)
      if (scell2[N_zone-1].star_num[j] == oldnum)
        {
        thisnum = j+1 ;
	if (t > t_start && t < t_end) fprintf(stdout, "Here? j=%4d\n", j) ;
        break ;
        }
   }

 previ_N_zone = N_zone ;

 intptr = candi_num ;

 if (thisnum == 0) thisnum = 1 ;  
 if (thisnum > scell2[N_zone-1].num_stars) 
     {
     thisnum = scell2[N_zone-1].num_stars ;  
     if (t > t_start && t < t_end) fprintf(stdout, "Here2? thisnum=%4d\n", thisnum) ;
     }
 thisnum1 = thisnum0 = thisnum ;

 if ( BD_dec >= 0.0) BD_dec += 6.0 ; /* 4.0 */
 if ( BD_dec <  0.0) BD_dec -= 6.0 ; /* 4.0 */

 i=1 ; 
 do{
  i_ra  = stars[scell2[N_zone-1].star_num[thisnum-i]-1].ra  ;
  i_dec = stars[scell2[N_zone-1].star_num[thisnum-i]-1].dec ;
  i_mag = stars[scell2[N_zone-1].star_num[thisnum-i]-1].mag ;
  
  ra_diff_new = fabs(*new_ra - i_ra) ;
  ra_diff_BD  = fabs(BD_ra  - i_ra) ;
  if (ra_diff_new > 180.0)  ra_diff_new = 360.0 - ra_diff_new ;
  if (ra_diff_BD  > 180.0)  ra_diff_BD  = 360.0 - ra_diff_BD  ;
  ra_diff_new *= cos(*new_dec*PI/180.0) * 3600.0 ;
  ra_diff_BD  *= cos(  BD_dec*PI/180.0) * 3600.0 ;

  if ( fabs(*new_dec) > 80.0 || fabs(BD_dec) > 80.0)   /* 88 */
     {
     ra_diff_BD  = 100.0 ;
     ra_diff_new = 100.0 ;
     }

  i_L[0] = cos(i_dec*PI/180.0)*cos(i_ra*PI/180.0) ;
  i_L[1] = cos(i_dec*PI/180.0)*sin(i_ra*PI/180.0) ;
  i_L[2] = sin(i_dec*PI/180.0) ;

  vec_dot = new->L[0]*i_L[0] + new->L[1]*i_L[1] + new->L[2]*i_L[2] ;
  if (vec_dot > 1.0)  vec_dot =  1.0 ;
  if (vec_dot < -1.0) vec_dot = -1.0 ;
  arc_dist = acos(vec_dot)*180.0/PI*3600.0 ;

  /* number = scell2[N_zone-1].star_num[thisnum+i] ; */

  if(fabs(arc_dist) < vici && fabs(i_mag - *new_mag) < mTOL) /* a possible match */
    {
    if (t > t_start && t < t_end) fprintf(stdout, "Here+-? thisnum=%4d  i=%3d\n", thisnum, i) ;
    thisnum1 = thisnum-i+1 ; 
    for(j=0;j<id_count;j++)
       if( id_star[j] == stars[scell2[N_zone-1].star_num[thisnum-i]-1].cat_num) 
         astaragain = 1 ;
    if (astaragain != 1)  /* no repeated star */    
       {
       *intptr  = stars[scell2[N_zone-1].star_num[thisnum-i]-1].cat_num  ;
       fprintf(outdist, "2  %15.6f %5d %15.8f \n", t, *intptr, arc_dist) ;
       (*count)++ ;
       intptr++ ;
       }
    }
  i++ ;
  astaragain = 0 ;
  if( (thisnum - i) < 0 )  /* break ;  */
	{
        thisnum = thisnum + scell2[N_zone-1].num_stars ;  
        if (t > t_start && t < t_end) fprintf(stdout, "Here3? thisnum=%4d\n", thisnum) ;
	}
  if (t > t_start && t < t_end) 
    {
    fprintf(stdout, "+++ ra_diff_BD = %9.2f, ra_diff_new = %9.2f\n",
           ra_diff_BD, ra_diff_new) ;
    fprintf(stdout, "  i = %3d  <<<<  i_con = %6.2f\n", 
	   i, scell2[N_zone-1].num_stars/2.0) ;
    }
  } while (  (  (fabs(ra_diff_BD)  < BD_limit*3600+1800.)    /* HERE */
         ||   (fabs(ra_diff_new) < FOV_limit*3600+1800.) ) 
         &&  i < (int) floor(scell2[N_zone-1].num_stars/2.0) ) ;

 thisnum = thisnum0 ;
 if (t > t_start && t < t_end) fprintf(stdout, "Here4? thisnum=%4d i=%3d\n", thisnum, i) ;
 if (thisnum == scell2[N_zone-1].num_stars)
           thisnum = 1 ; /* was 1 */  
 i=0 ;
 do{
  i_ra  = stars[scell2[N_zone-1].star_num[thisnum+i]-1].ra  ;
  i_dec = stars[scell2[N_zone-1].star_num[thisnum+i]-1].dec ;
  i_mag = stars[scell2[N_zone-1].star_num[thisnum+i]-1].mag ;

  if (t > t_start && t < t_end) 
    {
    fprintf(stdout, "2ND: t=%15.6f N_zone=%3d thisnum=%3d\n", t, N_zone, thisnum) ;
    fprintf(stdout, "i_ra = %9.2f, i_dec = %9.2f, i_mag = %6.2f\n",
	       i_ra, i_dec, i_mag) ;
    fprintf(stdout, "    star_num = %5d (%2d)\n", 
	   stars[scell2[N_zone-1].star_num[thisnum+i]-1].cat_num, i) ;
    }

  ra_diff_new = fabs(*new_ra - i_ra) ;
  ra_diff_new = fabs(*new_ra - i_ra)  ;
  ra_diff_BD  = fabs(BD_ra  - i_ra)  ;
  if (ra_diff_new > 180.0)  ra_diff_new = 360.0 - ra_diff_new ;
  if (ra_diff_BD  > 180.0)  ra_diff_BD  = 360.0 - ra_diff_BD  ;
  ra_diff_new *= cos(*new_dec*PI/180.0) * 3600.0 ;
  ra_diff_BD  *= cos(  BD_dec*PI/180.0) * 3600.0 ;

  if ( fabs(*new_dec) > 80.0 || fabs(BD_dec) > 80.0)   /* 88 */ 
     {
     ra_diff_BD  = 100.0 ;
     ra_diff_new = 100.0 ;
     }
  
  i_L[0] = cos(i_dec*PI/180.0)*cos(i_ra*PI/180.0) ;
  i_L[1] = cos(i_dec*PI/180.0)*sin(i_ra*PI/180.0) ;
  i_L[2] = sin(i_dec*PI/180.0) ;

  vec_dot = new->L[0]*i_L[0] + new->L[1]*i_L[1] + new->L[2]*i_L[2] ;
  if (vec_dot > 1.0)  vec_dot =  1.0 ;
  if (vec_dot < -1.0) vec_dot = -1.0 ;
  arc_dist = acos(vec_dot)*180.0/PI*3600.0 ;

  /* number = scell2[N_zone-1].star_num[thisnum-i] ; */

  if(fabs(arc_dist) < vici && fabs(i_mag - *new_mag) < mTOL) /* a possible match */
    {
    thisnum1 = thisnum+i+1 ;
    for(j=0;j<id_count;j++)
       if( id_star[j] == stars[scell2[N_zone-1].star_num[thisnum+i]-1].cat_num) 
         astaragain = 1 ;

    if( astaragain != 1)  /* no repeated star */    
      {
      *intptr  = stars[scell2[N_zone-1].star_num[thisnum+i]-1].cat_num  ;
      fprintf(outdist, "2  %15.6f %5d %15.8f \n", t, *intptr, arc_dist) ;
      (*count)++ ;
      intptr++ ;
      }
    }

  i++ ;
  astaragain = 0 ;
  if( (thisnum + i) >= scell2[N_zone-1].num_stars) /* break ; */
	{
        thisnum = thisnum - scell2[N_zone-1].num_stars ;  
        if (t > t_start && t < t_end) fprintf(stdout, "Here5? thisnum=%4d\n", thisnum) ;
	}

  } while (   ( (fabs(ra_diff_BD)  < BD_limit*3600+1800.)   /* HERE */
          ||  (fabs(ra_diff_new) < FOV_limit*3600+1800.) )  
          &&  i < (int) floor(scell2[N_zone-1].num_stars/2.0) ) ;

  thisnum = thisnum1 ;  

 return(1) ;
 }

int zone(dec)     /* find the zone where the BD is pointing */
double  dec ;
 {
    if (dec >= 80.0)                 return(1) ;
    if (dec <  80.0 && dec >=  75.0) return(2) ;
    if (dec <  75.0 && dec >=  70.0) return(3) ;
    if (dec <  70.0 && dec >=  65.0) return(4) ;
    if (dec <  65.0 && dec >=  60.0) return(5) ;
    if (dec <  60.0 && dec >=  55.0) return(6) ;
    if (dec <  55.0 && dec >=  50.0) return(7) ;
    if (dec <  50.0 && dec >=  45.0) return(8) ;
    if (dec <  45.0 && dec >=  40.0) return(9) ;
    if (dec <  40.0 && dec >=  35.0) return(10) ;
    if (dec <  35.0 && dec >=  30.0) return(11) ;
    if (dec <  30.0 && dec >=  25.0) return(12) ;
    if (dec <  25.0 && dec >=  20.0) return(13) ;
    if (dec <  20.0 && dec >=  15.0) return(14) ;
    if (dec <  15.0 && dec >=  10.0) return(15) ;
    if (dec <  10.0 && dec >=   5.0) return(16) ;
    if (dec <   5.0 && dec >=   0.0) return(17) ;
    if (dec <   0.0 && dec >=  -5.0) return(18) ;
    if (dec <  -5.0 && dec >= -10.0) return(19) ;
    if (dec < -10.0 && dec >= -15.0) return(20) ;
    if (dec < -15.0 && dec >= -20.0) return(21) ;
    if (dec < -20.0 && dec >= -25.0) return(22) ;
    if (dec < -25.0 && dec >= -30.0) return(23) ;
    if (dec < -30.0 && dec >= -35.0) return(24) ;
    if (dec < -35.0 && dec >= -40.0) return(25) ;
    if (dec < -40.0 && dec >= -45.0) return(26) ;
    if (dec < -45.0 && dec >= -50.0) return(27) ;
    if (dec < -50.0 && dec >= -55.0) return(28) ;
    if (dec < -55.0 && dec >= -60.0) return(29) ;
    if (dec < -60.0 && dec >= -65.0) return(30) ;
    if (dec < -65.0 && dec >= -70.0) return(31) ;
    if (dec < -70.0 && dec >= -75.0) return(32) ;
    if (dec < -75.0 && dec >= -80.0) return(33) ;
    if (dec < -80.0 )                return(34) ;

 return(-1)  ;
 }

int  multi_obs(t, new, several_num, count)
struct VecStar *new ;
int    *several_num, count ;
double t ;   /* t : not used */
 {
 int    i, counted=0, *intptr, iptr, *cptr ; 
 int    dis_compare() ;
 double mag_diff ;

 intptr = (int *) malloc(count*sizeof(int)) ;  /* Ns */
 if ( intptr == NULL) {
   fprintf(stdout, " FAIL : malloc intptr (starID_dm::multi_obs) \n") ;
   exit(5) ; }

 cptr = intptr ;
 for(i=0;i<count;i++)
   {
   mag_diff=(new->mag-stars[*(several_num+i)-1].mag) ;
   if (fabs(mag_diff) < mTOL)
      {
      counted++ ;
      *cptr = stars[*(several_num+i)-1].cat_num ;
      cptr++ ;
      }
   }

 switch(counted)
   {
   case 0 :  iptr = dis_compare(t, new, several_num, count) ;
    	     return(iptr) ;
   case 1 :  iptr = (*cptr-1) ;
  	     return(iptr) ;
   default:  iptr = dis_compare(t, new, intptr, counted) ;
	     return(iptr) ;
   }

 fprintf(stdout,"starID_dm():find():multi_obs() didn't finish normally.\n") ;
 exit(6) ;
 }

int dis_compare(t, new, intptr, counted)
struct VecStar *new ;
int    counted, *intptr ;
double t ;  /* t : not used */
 {
 int     i, smallest ;
 double  i_L[3], vec_dot, arc_dist, arc_dist_s ;

 i_L[0] = cos(stars[*intptr-1].dec*PI/180.0)*cos(stars[*intptr-1].ra*PI/180.0) ;
 i_L[1] = cos(stars[*intptr-1].dec*PI/180.0)*sin(stars[*intptr-1].ra*PI/180.0) ;
 i_L[2] = sin(stars[*intptr-1].dec*PI/180.0) ;
 vec_dot = new->L[0]*i_L[0] + new->L[1]*i_L[1] + new->L[2]*i_L[2] ;
 if (vec_dot > 1.0)  vec_dot =  1.0 ;
 if (vec_dot < -1.0) vec_dot = -1.0 ;
 arc_dist_s = acos(vec_dot) ;
 smallest = *intptr ;
     
 for(i=1;i<counted;i++)
   {
   i_L[0] = cos(stars[*(intptr+i)-1].dec*PI/180.0)*cos(stars[*(intptr+i)-1].ra*PI/180.0) ;
   i_L[1] = cos(stars[*(intptr+i)-1].dec*PI/180.0)*sin(stars[*(intptr+i)-1].ra*PI/180.0) ;
   i_L[2] = sin(stars[*(intptr+i)-1].dec*PI/180.0) ;

   vec_dot = new->L[0]*i_L[0] + new->L[1]*i_L[1] + new->L[2]*i_L[2] ;
   if (vec_dot > 1.0)  vec_dot =  1.0 ;
   if (vec_dot < -1.0) vec_dot = -1.0 ;
   arc_dist = acos(vec_dot) ;

   if( fabs(arc_dist) < fabs(arc_dist_s))
      {
      smallest = *intptr+i ;
      arc_dist_s = arc_dist ;
      }
   }
 return(smallest) ; 
}

