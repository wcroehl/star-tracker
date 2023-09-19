#define NUMs 6
starID_pm(t, q0, nstar, id_count, bore_body, 
      outi, outr, outo, outmag, ixy, x, y, xmag, outes1, outes2)
FILE    *outi, *outr, *outo, *outmag, *ixy, *outes1, *outes2; /*added outes1 and outes2, 7/28/23*/
int     *id_count, nstar ;  /* the number of measured stars */
double  t, q0[4], bore_body[3], x[NUMs], y[NUMs], xmag[NUMs] ;
{
  short int  virtual_star() ;
  int     i, j, k, k2, ij, ii, jj, i0, j0, esignal, rtn,
          nro, nb_cell, cellnum, starnum,   
          cell(), sort() ;  
  double  al_rad, de_rad,   /* RA and Dec in radian */
          est_alp, est_del, norm, bore_iner[3], pos[3], Ap[3][3],
          BD[3], L[3], ang_dist ;   
  void    cartesian(), obs_ang(),  rearrange(),
          on6stars(),  on5stars(), on4stars(), on3stars(), 
          compare0(),  compare1() ;

  ij = 0 ;
  rtn = 0 ;
  *id_count = 0 ;

  q_to_A(q0, Ap) ;

  for(i=0;i<3;i++)   bore_iner[i] = 0.0 ;
  for(i=0;i<3;i++) for(j=0;j<3;j++)
        bore_iner[i] += Ap[j][i] * bore_body[j] ;

  est_del = asin(bore_iner[2])/PI*180.0 ;              
  est_alp = atan2(bore_iner[1], bore_iner[0])/PI*180.0; /* [][1],[][0] */
  if(bore_iner[1] < 0.0) est_alp += 360.0 ;

  nb_cell = cell(est_alp, est_del)-1 ; /* cell num. associated with BD */
                               /* -1 added since array start with zero */

  if (nb_cell == -1)
     {
     fprintf(stderr,"starID_pm():cell() didn't finish normally.\n") ;
     return (6) ;
     }

  cartesian(est_alp*PI/180, est_del*PI/180, BD) ;
  for (i=0 ; i < adjcell[nb_cell].num_adj ; i++)
      {
      cellnum = adjcell[nb_cell].cell_num[i]-1 ; /* read adjacent cell # */ 
      for(j=0 ; j < scell[cellnum].num_stars ; j++)
          {                 /* read all stars in the cell */ 
          starnum = scell[cellnum].star_num[j]-1 ;  /* catalog # */
          al_rad = stars[starnum].ra  / 180 * PI ;   /* RA  in radian */
          de_rad = stars[starnum].dec / 180 * PI ;   /* Dec in radian */
          cartesian(al_rad,de_rad,L) ;  /* get star pos. in cartesian */
 
          BLI[ij].IBL = BD[0]*L[0] + BD[1]*L[1] + BD[2]*L[2] ;
                 /* dot product of BD and L --> the angle of cosine */
          BLI[ij].starnum = starnum+1 ;  /* add 1 from array number */
          for (k=0;k<3;k++) BLI[ij].L[k] = L[k] ;
          BLI[ij].mag  = stars[starnum].mag ; 
          ij++ ; /* the total # of candidate stars from star catalog */
          }  /* end of individual star input */
      }  /*  end of <for(for)> statement */

   rtn = sort(ij) ;  /* structure BLI reorganized with the angular
                   distance between BD and star pos. smaller first */
   if (rtn != 1)
   {
   fprintf(stderr,"starID_pm():sort() didn't finish normally.\n") ;
   exit(6) ;
   }

   nro = nstar - 1 ; /* the possible # of rotations of measured stars */

   k = 0 ; esignal = 3 ; k2 = 0 ; /* esignal = 3 : not IDed */

   while (esignal == 3 && k < nro) {
      for(i=0;i<NUMs;i++) 
         {
         NEWBLI[i].mag = 0.0 ;
         NEWBLI[i].id = -999 ;
         NEWBLI[i].num = 0 ;  
         }
      compare0(t, 0, 1, &ii, &jj, ij, &esignal) ;

      fprintf(outes1, "%15.6f %15.6f %15.6f %15.6f \n", t, ii, jj, esignal);
      /*added the above fprintf function, 7/28/23*/

      i0 = ii ;
      j0 = jj ;

      if (esignal == 3) { obs_ang(nstar, x, y, xmag) ; k++ ; }
      } 

   if (esignal < 3) 
      do {
         switch(nstar)
            { 
            case 6 : on6stars(t, &jj,ij,&esignal) ;
                     break ;
            case 5 : on5stars(t, &jj,ij,&esignal) ;
                     break ;
            case 4 : on4stars(t, &jj,ij,&esignal) ;

                     fprintf(outes2, "%15.6f %15.6f %15.6f \n", t, jj, esignal);
                     /*added the above fprintf function, 7/28/23*/
                     break ;
            case 3 : on3stars(t, &jj,ij,&esignal) ;
            }
         if (esignal == 1)  { ii = j0 ; jj = i0 ; }  
         while (esignal == 3 & k < nro && i0 < ij && j0 < ij) 
            {
            ii = i0 ;  jj = j0 ; 
            for(i=0;i<NUMs;i++) 
               {
               NEWBLI[i].mag = 0.0 ;
               NEWBLI[i].id = -999 ;
               NEWBLI[i].num = 0 ;
               }
            compare1(t, 0, 1, &ii, &jj, ij, &esignal) ; 
            while (esignal == 3 && k < nro && ii == ij && jj == ij)
               {
               obs_ang(nstar, x, y, xmag) ; k++ ;
               compare0(t, 0, 1, &ii, &jj, ij, &esignal) ; 
               } 
            } 
         /*  More search for base stars when nstar = 6 */
         if (esignal == 3 && k == nro && nstar == 6 && k2 == 0 
                    && ii == ij && jj == ij) 
               {
               rearrange(nstar, x, y, xmag) ; k2=2 ; k=1 ;
               }
         while (esignal == 3 && k < nro && nstar == 6 && k2 > 1 
                    && ii == ij && jj == ij)
               { 
               obs_ang(nstar, x, y, xmag) ; k++ ;
               compare0(t, 0, 1, &ii, &jj, ij, &esignal) ; 
               k2 = 5 ;
               }
         if (k2 != 0) k2 = -1 ; 
         
         i0 = ii ;  j0 = jj ;
    } while (esignal < 3) ;  

if (esignal == 3) nstar = 0 ;

for(i=0;i<nstar;i++) 
   if (NEWBLI[i].mag != 0.0) 
     {
       crf[(*id_count)]->mag   = NEWBLI[i].mag ;
       for(j=0;j<3;j++) pos[j] = NEWBLI[i].L[j] ;
       crf[(*id_count)]->num   = NEWBLI[i].num ;
       crf[(*id_count)]->id    = NEWBLI[i].id ;

       norm = NEWBLI[i].L[0]*pos[0]+NEWBLI[i].L[1]*pos[1]+NEWBLI[i].L[2]*pos[2] ;

       if (norm >  1.0)  norm =  1.0 ;
       if (norm < -1.0)  norm = -1.0 ;
       ang_dist = acos(norm)*180.0/PI*3600.0; 
	      /* angular dist. between before&after aberration correction */

       /* aberration not simulated */
       for (j=0;j<3;j++) crf[(*id_count)]->L[j] = pos[j] ; 
					 /* for the aberration correction */
       for(j=0;j<3;j++) crf[(*id_count)]->L[j] = NEWBLI[i].L[j] ;

       for(j=0;j<3;j++) m_star[(*id_count)]->L[j] = m_star[i]->L[j] ;
       x[(*id_count)] = x[i] ;
       y[(*id_count)] = y[i] ;

       /* ccd frames to Body-fixed frame  */
       for(j=0;j<3;j++) b_star[(*id_count)]->L[j] = 0.0 ;
       for(j=0;j<3;j++) for(k=0;k<3;k++)
          b_star[(*id_count)]->L[j] += T_B[j][k] * m_star[(*id_count)]->L[k] ;
       
       norm = sqrt(b_star[(*id_count)]->L[0]*b_star[(*id_count)]->L[0]
            +      b_star[(*id_count)]->L[1]*b_star[(*id_count)]->L[1]
            +      b_star[(*id_count)]->L[2]*b_star[(*id_count)]->L[2]) ;

       for(j=0;j<3;j++) b_star[(*id_count)]->L[j] /= norm ; 

       b_star[(*id_count)]->mag  = m_star[(*id_count)]->mag  = m_star[i]->mag ; 
       xmag[(*id_count)] = xmag[i] ;

       fprintf(outmag," %15.3f  1  %6d  %6.2f %6.2f %8.4f\n",
          t, crf[(*id_count)]->num, crf[(*id_count)]->mag, b_star[(*id_count)]->mag, ang_dist) ;

       fprintf(outi,"%15.6f  %20.15f %20.15f %20.15f %10.2f       %10d\n",
           t, crf[(*id_count)]->L[0], crf[(*id_count)]->L[1], crf[(*id_count)]->L[2],
           crf[(*id_count)]->mag, crf[(*id_count)]->num) ;   

       fprintf(outo,"%15.6f  %20.15f %20.15f %20.15f %10.2f\n",
           t, m_star[i]->L[0], m_star[i]->L[1], m_star[i]->L[2],
           m_star[i]->mag) ;   
       
       fprintf(ixy, "%15.6f  %10.2f %10.2f %10.2f %10.2f %10d\n",
           t, x[i], y[i], xmag[i], m_star[i]->mag, crf[(*id_count)]->num) ;   

       fprintf(outr,"%15.6f  %20.15f %20.15f %20.15f %10.2f\n",
           t, b_star[(*id_count)]->L[0], b_star[(*id_count)]->L[1], b_star[(*id_count)]->L[2],
           b_star[(*id_count)]->mag) ;   
       (*id_count)++ ;
     } /* the end of if (NEWBLI[i].mag != 0.0) */
return(1) ;
}       

/*  find the primary cell number (from the boresight vector) 
    corresponding to the estimated RA and Dec   */
cell(RA, dec)
double  RA, dec ;
{
    int       n, j, No, cellsu, jmax ;
    double    dec_dif, RA_dif, de ;

    No = 14 ; /* 15 declinamtion zones */

    if (dec <  90.00 && dec >=  83.79)  cellsu =  1 ;
    if (dec <  83.79 && dec >=  71.38)  cellsu =  5 ;
    if (dec <  71.38 && dec >=  58.97)  cellsu =  9 ;
    if (dec <  58.97 && dec >=  46.55)  cellsu = 13 ;
    if (dec <  46.55 && dec >=  34.14)  cellsu = 17 ;
    if (dec <  34.14 && dec >=  21.72)  cellsu = 21 ;
    if (dec <  21.72 && dec >=   9.31)  cellsu = 25 ;
    if (dec <   9.31 && dec >=  -3.10)  cellsu = 29 ;
    if (dec <  -3.10 && dec >= -15.52)  cellsu = 27 ;
    if (dec < -15.52 && dec >= -27.93)  cellsu = 23 ;
    if (dec < -27.93 && dec >= -40.34)  cellsu = 19 ;
    if (dec < -40.34 && dec >= -52.76)  cellsu = 15 ;
    if (dec < -52.76 && dec >= -65.17)  cellsu = 11 ;
    if (dec < -65.17 && dec >= -77.59)  cellsu =  7 ;
    if (dec < -77.59 && dec >= -90.00)  cellsu =  3 ;

    dec_dif = 12.4136 ;  /* 12.4136 = (83.79+90.00)/14 */
    RA_dif = 360./cellsu ;
    de = 90. - dec ;

    j = (int) ceil(RA/RA_dif - 0.5) ;
    jmax = (int) ceil(360./RA_dif - 0.5) ;

    if (de < 93.1001295051)
        n = (int) 2*ceil(de/dec_dif - 0.5) ;
    else
        n = 2*No + 1 - (int) 2*ceil(de/dec_dif - 0.5) ;
    if ( (n*n+j+1) == (n*n+jmax+1) )  return(n*n+1) ;
    return(n*n+j+1) ; 
}

sort(ij)  /* sorting BLI from smallest ang. dist. */
int ij ;
{
    int i, j ;
    struct list temporary ;

    for(i=0;i<ij;i++) { 
       for(j=i+1;j<ij;j++)  /* smaller BLI first */
          if(BLI[j].IBL > BLI[i].IBL) 
            { 
            temporary.starnum = BLI[i].starnum ;
            temporary.IBL     = BLI[i].IBL     ;
            temporary.L[0]    = BLI[i].L[0]    ;
            temporary.L[1]    = BLI[i].L[1]    ;
            temporary.L[2]    = BLI[i].L[2]    ;
            temporary.mag     = BLI[i].mag     ;

            BLI[i].starnum = BLI[j].starnum ;
            BLI[i].IBL     = BLI[j].IBL     ;
            BLI[i].L[0]    = BLI[j].L[0]    ;
            BLI[i].L[1]    = BLI[j].L[1]    ;
            BLI[i].L[2]    = BLI[j].L[2]    ;
            BLI[i].mag     = BLI[j].mag     ;

            BLI[j].starnum = temporary.starnum ;
            BLI[j].IBL     = temporary.IBL     ;
            BLI[j].L[0]    = temporary.L[0]    ;
            BLI[j].L[1]    = temporary.L[1]    ;
            BLI[j].L[2]    = temporary.L[2]    ;
            BLI[j].mag     = temporary.mag     ;
            } 
       } 
return(1) ;
} 

void on6stars(t, jj, ij, esignal)
int     *jj, ij, *esignal ;
double  t ;
{
   int    k, sign ;
   void   compare2(), compare3(), compare4() ;
   double tolerance ;

   tolerance = TOL6 ; 
   compare2(0, 2, &(*jj), ij, &sign, tolerance) ;  /* 3rd star search */
   if(sign != 2) compare3(1, 2, *jj, &sign, tolerance) ;
   while(sign == 2  && *jj < ij)    
      {
      compare4(0, 2, &(*jj), ij, &sign, tolerance) ;
      if(sign != 2) compare3(1, 2, *jj, &sign, tolerance) ;
      }
   if(sign != 2)   /* 3rd star matched, 4th star search */
      {
      compare2(0, 3, &(*jj), ij, &sign, tolerance) ;
      if(sign != 2) compare3(1, 3, *jj, &sign, tolerance) ;
      while(sign == 2  && *jj < ij) 
         {
         compare4(0, 3, &(*jj), ij, &sign, tolerance) ;
         if(sign != 2) compare3(1, 3, *jj, &sign, tolerance);
         }
      if(sign != 2)  /* 3rd, 4th matched, 5th star search */ 
         { 
         compare2(0, 4, &(*jj), ij, &sign, tolerance);
         if(sign != 2) compare3(1, 4, *jj, &sign, tolerance) ;
         while(sign == 2  && *jj < ij) 
            {
            compare4(0, 4, &(*jj), ij, &sign, tolerance) ;
            if(sign != 2) compare3(1, 4, *jj, &sign, tolerance);
            }
         if(sign != 2) /* 3rd, 4th, 5th matched, 6th search */ 
            { 
            compare2(0, 5, &(*jj), ij, &sign, tolerance);
            if(sign != 2) compare3(1, 5, *jj, &sign, tolerance) ;
            while(sign == 2  && *jj < ij) 
               {
               compare4(0, 5, &(*jj), ij, &sign, tolerance) ;
               if(sign != 2) compare3(1, 5, *jj, &sign, tolerance);
               }
            if(sign != 2) *esignal = 6 ;   /* 6 stars */
            if(sign == 2) *esignal = 5 ;   /* 5 stars (0,1,2,3,4) */ 
            return ;
            }
         if(sign == 2) /* 3rd, 4th matched, 5th not, 6th search */ 
            { 
            compare2(0, 5, &(*jj), ij, &sign, tolerance);
            if(sign != 2) compare3(1, 5, *jj, &sign, tolerance) ;
            while(sign == 2  && *jj < ij) 
               {
               compare4(0, 5, &(*jj), ij, &sign, tolerance) ;
               if(sign != 2) compare3(1, 5, *jj, &sign, tolerance);
               }
            if(sign != 2)  /* 3rd, 4th, 6th matched, 5th not */
               { 
               *esignal = 5 ;         /* 5 stars (0,1,2,3,5) */
               return ; 
               }
            if(sign == 2)  /* 3rd, 4th matched, 5th, 6th not */
               { 
               *esignal = 4 ;         /* 4 stars (0,1,2,3)   */
               return ; 
               }
            }
         } 
      if(sign == 2)  /* 3rd matched, 4th not, 5th star search */ 
         {
         compare2(0, 4, &(*jj), ij, &sign, tolerance);
         if(sign != 2) compare3(1, 4, *jj, &sign, tolerance) ;
         while(sign == 2  && *jj < ij) 
            {
            compare4(0, 4, &(*jj), ij, &sign, tolerance) ;
            if(sign != 2) compare3(1, 4, *jj, &sign, tolerance);
            }
         if(sign != 2) /* 3rd OK, 4th not, 5th matched, 6th search */ 
            { 
            compare2(0, 5, &(*jj), ij, &sign, tolerance);
            if(sign != 2) compare3(1, 5, *jj, &sign, tolerance) ;
            while(sign == 2  && *jj < ij) 
               {
               compare4(0, 5, &(*jj), ij, &sign, tolerance) ;
               if(sign != 2) compare3(1, 5, *jj, &sign, tolerance);
               }
            if(sign != 2)    *esignal = 5 ;  /* 5 stars (0,1,2,4,5) */
            if(sign == 2)    *esignal = 4 ;  /* 4 stars (0,1,2,4)   */ 
            return ;
            }
         if(sign == 2) /* 3rd OK, 4th not, 5th not, 6th search */ 
            { 
            compare2(0, 5, &(*jj), ij, &sign, tolerance);
            if(sign != 2) compare3(1, 5, *jj, &sign, tolerance) ;
            while(sign == 2  && *jj < ij) 
               {
               compare4(0, 5, &(*jj), ij, &sign, tolerance) ;
               if(sign != 2) compare3(1, 5, *jj, &sign, tolerance);
               }
            if(sign != 2)  /* 3rd OK, 4th not, 5th not, 6th OK */
               { 
               *esignal = 4 ;       /* 4 stars (0,1,2,5) */
               return ;
               }
            if(sign == 2)  /* 3rd OK, 4th, 5th, 6th not matched */ 
               {
               /* 3 stars among 6 or (undetected) more stars in the FOV
                  may easily cause ID error */
                  sign = 3 ; 
               }
            }
         }
      } 

   if(sign == 2)   /* 3rd not matched, 4th star search */
      {
      compare2(0, 3, &(*jj), ij, &sign, tolerance) ;
      if(sign != 2) compare3(1, 3, *jj, &sign, tolerance) ;
      while(sign == 2  && *jj < ij) 
         {
         compare4(0, 3, &(*jj), ij, &sign, tolerance) ;
         if(sign != 2) compare3(1, 3, *jj, &sign, tolerance);
         }
      if(sign != 2)  /* 3rd not matched, 4th OK, 5th search */ 
         {
         compare2(0, 4, &(*jj), ij, &sign, tolerance);
         if(sign != 2) compare3(1, 4, *jj, &sign, tolerance) ;
         while(sign == 2  && *jj < ij) 
            {
            compare4(0, 4, &(*jj), ij, &sign, tolerance) ;
            if(sign != 2) compare3(1, 4, *jj, &sign, tolerance);
            }
         if(sign != 2) /* 3rd not, 4th, 5th matched, 6th search */ 
            {
            compare2(0, 5, &(*jj), ij, &sign, tolerance);
            if(sign != 2) compare3(1, 5, *jj, &sign, tolerance);
            while(sign == 2  && *jj < ij) 
                {
                compare4(0, 5, &(*jj), ij, &sign, tolerance) ;
                if(sign != 2) compare3(1, 5, *jj, &sign, tolerance) ;
                }
            if (sign != 2)   *esignal = 5 ; /* 5 stars (0,1,3,4,5) */
            if (sign == 2)   *esignal = 4 ; /* 4 stars (0,1,3,4)   */
            return ; 
            }
         if(sign == 2) /* 3rd not, 4th OK, 5th not, 6th search */ 
            {
            compare2(0, 5, &(*jj), ij, &sign, tolerance);
            if(sign != 2) compare3(1, 5, *jj, &sign, tolerance);
            while(sign == 2  && *jj < ij) 
                {
                compare4(0, 5, &(*jj), ij, &sign, tolerance) ;
                if(sign != 2) compare3(1, 5, *jj, &sign, tolerance) ;
                }
            if (sign != 2) /* 3rd not, 4th OK, 5th not, 6th OK */ 
                {
                *esignal = 4 ;    /* 4 stars (0,1,3,5) */
                return ; 
                }
            if (sign == 2) sign = 3 ; /* 3rd not, 4th OK, 5th, 6th not */ 
            }
         }
      if(sign == 2)  /* 3rd not matched, 4th not, 5th search */ 
         {
         compare2(0, 4, &(*jj), ij, &sign, tolerance);
         if(sign != 2) compare3(1, 4, *jj, &sign, tolerance) ;
         while(sign == 2  && *jj < ij) 
            {
            compare4(0, 4, &(*jj), ij, &sign, tolerance) ;
            if(sign != 2) compare3(1, 4, *jj, &sign, tolerance);
            }
         if(sign != 2) /* 3rd, 4th not, 5th matched, 6th search */ 
            {
            compare2(0, 5, &(*jj), ij, &sign, tolerance);
            if(sign != 2) compare3(1, 5, *jj, &sign, tolerance);
            while(sign == 2  && *jj < ij) 
                {
                compare4(0, 5, &(*jj), ij, &sign, tolerance) ;
                if(sign != 2) compare3(1, 5, *jj, &sign, tolerance) ;
                }
            if (sign != 2)  /* 3rd, 4th not, 5th, 6th matched */
                {
                *esignal = 4 ;        /* 4 stars (0,1,4,5) */
                return ;
                }
            if (sign == 2)  sign = 3 ; /* 3rd, 4th not, 5th matched, 6th not */ 
            }
         if(sign == 2) /* 3rd, 4th, 5th not matched, how about 6th? */ 
            {
            compare2(0, 5, &(*jj), ij, &sign, tolerance);
            if(sign != 2) compare3(1, 5, *jj, &sign, tolerance) ;
            while(sign == 2  && *jj < ij) 
               {
               compare4(0, 5, &(*jj), ij, &sign, tolerance) ;
               if(sign != 2) compare3(1, 5, *jj, &sign, tolerance);
               }
            if(sign != 2) sign = 3 ; /* 3rd, 4th, 5th not matched, 6th OK */ 
            /* if(sign == 2) ; */
            } 
         }
      }         
   if (*esignal == 0)
      {
      for(k=0;k<3;k++)
         NBLI.L[k] = NEWBLI[0].L[k] ;
         NBLI.num = NEWBLI[0].num ;
         NBLI.mag   = NEWBLI[0].mag ; 
         NBLI.id = NEWBLI[0].id ;

      for(k=0;k<3;k++)
         NEWBLI[0].L[k] = NEWBLI[1].L[k] ;
         NEWBLI[0].num = NEWBLI[1].num ;
         NEWBLI[0].mag   = NEWBLI[1].mag ;
         NEWBLI[0].id = NEWBLI[1].id ;

      for(k=0;k<3;k++)
         NEWBLI[1].L[k] = NBLI.L[k] ;
         NEWBLI[1].num = NBLI.num ;
         NEWBLI[1].mag   = NBLI.mag ;
         NEWBLI[1].id = NBLI.id ;

      (*esignal)++ ;
      return ;
      }  
   if (*esignal <= 2)  
      {  
      *esignal = 3 ; 
      return ;
      }
return ;
}

void on5stars(t, jj, ij, esignal)
int    *jj, ij, *esignal ;
double t ;
{
   int    k, sign ; 
   void   compare2(), compare3(), compare4() ;
   double tolerance ;

   tolerance = TOL6 ;
   compare2(0, 2, &(*jj), ij, &sign, tolerance) ;  /* 3rd star searching */ 
   if(sign != 2) compare3(1, 2, *jj, &sign, tolerance) ;
   while(sign == 2  && *jj < ij) 
      {
      compare4(0, 2, &(*jj), ij, &sign, tolerance) ;
      if(sign != 2) compare3(1, 2, *jj, &sign, tolerance) ;
      }

   if(sign != 2)          /* 3rd star matched, 4th star search */ 
      {
      compare2(0, 3, &(*jj), ij, &sign, tolerance) ;
      if(sign != 2) compare3(1, 3, *jj, &sign, tolerance) ;
      while(sign == 2  && *jj < ij) 
          {
          compare4(0, 3, &(*jj), ij, &sign, tolerance) ;
          if(sign != 2) compare3(1, 3, *jj, &sign, tolerance);
          }

      if(sign != 2)       /* 3rd, 4th matched, 5th star search */ 
          { 
          compare2(0, 4, &(*jj), ij, &sign, tolerance);
          if(sign != 2) compare3(1, 4, *jj, &sign, tolerance) ;
          while(sign == 2  && *jj < ij) 
              {
	      compare4(0, 4, &(*jj), ij, &sign, tolerance) ;
	      if(sign != 2) compare3(1, 4, *jj, &sign, tolerance);
	      }
 	  if(sign != 2)   *esignal = 5 ;  /* 5 stars (all)     */
	  if(sign == 2)   *esignal = 4 ;  /* 4 stars (0,1,2,3) */ 
          return ;
          }
      if(sign == 2)   /* 3rd matched, 4th not, 5th star search */  
          { 
          compare2(0, 4, &(*jj), ij, &sign, tolerance);
          if(sign != 2) compare3(1, 4, *jj, &sign, tolerance) ;
          while(sign == 2  && *jj < ij) 
              {
	      compare4(0, 4, &(*jj), ij, &sign, tolerance) ;
              if(sign != 2) compare3(1, 4, *jj, &sign, tolerance);
	      }
          if(sign != 2)  /* 3rd OK, 4th not, 5th matched */ 
              { 
              *esignal = 4 ;     /* 4 stars (0,1,2,4) */
              return ; 
              }
          if(sign == 2)  /* 3rd OK, 4th, 5th not matched */ 
              { 
              sign = 3 ;
              }
          }
      }

   if(sign == 2)   /* 3rd star not matched, 4th star search */
      { 
      compare2(0, 3, &(*jj), ij, &sign, tolerance) ;
      if(sign != 2) compare3(1, 3, *jj, &sign, tolerance);
      while(sign == 2  && *jj < ij) 
          {
          compare4(0, 3, &(*jj), ij, &sign, tolerance) ;
          if(sign != 2) compare3(1, 3, *jj, &sign, tolerance);
          }
      if(sign != 2)   /* 3rd not matched, 4th matched, 5th search */ 
          { 
          compare2(0, 4, &(*jj), ij, &sign, tolerance);
          if(sign != 2) compare3(1, 4, *jj, &sign, tolerance);
          while(sign == 2  && *jj < ij) 
             {
             compare4(0, 4, &(*jj), ij, &sign, tolerance) ;
             if(sign != 2) compare3(1, 4, *jj, &sign, tolerance) ;
             }
	  if (sign != 2)  /* 3rd not matched, 4th 5th matched */ 
             { 
             *esignal = 4 ;    /* 4 stars (0,1,3,4) */ 
             return ; 
             }
	  if (sign == 2)  /* 3rd not, 4th OK, 5th not  */
             {
             sign = 3 ;
	     }
          }
      if(sign == 2)  /* 3rd 4th not matched. how about 5th  */ 
            {
            compare2(0, 4, &(*jj), ij, &sign, tolerance);
            if(sign != 2) compare3(1, 4, *jj, &sign, tolerance) ;
            while(sign == 2  && *jj < ij) 
               {
               compare4(0, 4, &(*jj), ij, &sign, tolerance) ;
               if(sign != 2) compare3(1, 4, *jj, &sign, tolerance);
               }
            if(sign != 2) /* 3rd, 4th not matched, 5th matched */ 
               {
               sign = 3 ;
               } 
            }
   }

   if (*esignal == 0) 
      {
      for(k=0;k<3;k++) 
         NBLI.L[k] = NEWBLI[0].L[k] ;
         NBLI.num = NEWBLI[0].num ;
         NBLI.mag   = NEWBLI[0].mag ;
         NBLI.id = NEWBLI[0].id ;

      for(k=0;k<3;k++) 
         NEWBLI[0].L[k] = NEWBLI[1].L[k] ;
         NEWBLI[0].num = NEWBLI[1].num ;
         NEWBLI[0].mag   = NEWBLI[1].mag ;
         NEWBLI[0].id = NEWBLI[1].id ;

      for(k=0;k<3;k++) 
         NEWBLI[1].L[k] = NBLI.L[k] ;
         NEWBLI[1].num = NBLI.num ;
         NEWBLI[1].mag   = NBLI.mag ;
         NEWBLI[1].id = NBLI.id ;

      (*esignal)++ ;
      return ;
      } ;  
   if (*esignal <= 2)  
      { 
      *esignal = 3 ; 
      return ;
      }
return ;
}

void on4stars(t, jj, ij, esignal)
int    *jj, *esignal, ij ;
double t ;
{
   int    k, sign ; 
   void   compare2(), compare3(), compare4() ;
   double tolerance ;
  
   tolerance = TOL3 ;
   compare2(0, 2, &(*jj), ij, &sign, tolerance) ; /* 3rd star search */
   if(sign != 2) compare3(1, 2, *jj, &sign, tolerance) ;
   while(sign == 2  && *jj < ij) 
      {
      compare4(0, 2, &(*jj), ij, &sign, tolerance) ;
      if(sign != 2) compare3(1, 2, *jj, &sign, tolerance) ;
      }

   if(sign != 2)  /* 3rd star matched, 4th star search */ 
      {
      compare2(0, 3, &(*jj), ij, &sign, tolerance) ;
      if(sign != 2) compare3(1, 3, *jj, &sign, tolerance) ;
      while(sign == 2  && *jj < ij) 
         {
         compare4(0, 3, &(*jj), ij, &sign, tolerance) ;
	 if(sign != 2) compare3(1, 3, *jj, &sign, tolerance);
	 }
      if(sign != 2) /* 3rd, 4th star matched */
         { 
         *esignal = 4 ;   /* 4 stars (all) */
         return ; 
         }
      if(sign == 2) /* 3rd OK, 4th not */
         {
         sign = 3 ;
         }
      }
   if(sign == 2)  /* 3rd star not matched, 4th star search */ 
      {
      compare2(0, 3, &(*jj), ij, &sign, tolerance) ;
      if(sign != 2) compare3(1, 3, *jj, &sign, tolerance) ;
      while(sign == 2  && *jj < ij) 
         {
         compare4(0, 3, &(*jj), ij, &sign, tolerance) ;
	 if(sign != 2) compare3(1, 3, *jj, &sign, tolerance);
	 }
      if(sign != 2) /* 3rd not match, 4th star matched */
         {
         sign = 3 ;
         }
      }
         
   if (*esignal == 0) 
      {
      for(k=0;k<3;k++) 
         NBLI.L[k] = NEWBLI[0].L[k] ;
      NBLI.num = NEWBLI[0].num ;
      NBLI.mag   = NEWBLI[0].mag ;
      NBLI.id = NEWBLI[0].id ;
   
      for(k=0;k<3;k++) 
         NEWBLI[0].L[k] = NEWBLI[1].L[k] ;
      NEWBLI[0].num = NEWBLI[1].num ;
      NEWBLI[0].mag   = NEWBLI[1].mag ;
      NEWBLI[0].id = NEWBLI[1].id ;
 
      for(k=0;k<3;k++) 
         NEWBLI[1].L[k] = NBLI.L[k] ;
      NEWBLI[1].num = NBLI.num ;
      NEWBLI[1].mag   = NBLI.mag ;
      NEWBLI[1].id = NBLI.id ;
   
      (*esignal)++ ;
      return ;
      } ;   
   if (*esignal <= 2)  
      { 
      *esignal = 3 ; 
      return ;
      }
return ;
}

void on3stars(t, jj, ij, esignal)
int    *jj, *esignal, ij ;
double t ; 
{
   int    k, sign, ambi, sym_test() ;
   void   compare2(), compare3(), compare4() ;
   double tolerance, dist_tol ;

   tolerance = TOL3 ;
   dist_tol  = DIST_TOL ;
   compare2(0, 2, &(*jj), ij, &sign, tolerance) ; /* 3rd star search */
   if(sign != 2) compare3(1, 2, *jj, &sign, tolerance) ;
   while(sign == 2  && *jj < ij)
      {
      compare4(0, 2, &(*jj), ij, &sign, tolerance) ;
      if(sign != 2) compare3(1, 2, *jj, &sign, tolerance) ;
      }
   if(sign != 2)
      {
      ambi = sym_test(dist_tol) ; 
      if (ambi == 0) { *esignal = 7 ;  /* 3 stars (all) */ } 
      else             *esignal = 3 ;
      }
   if(sign == 2) /* 3rd not matched. */
      {
      *esignal = 3 ;
      }
return ; 
}

/* compare the (cosine) of the interangles of pairs
   of BLI(catalog) pairs and meas_ang */
void compare0(t, n, m, ii, jj, ij, esignal)
int     n, m, *ii, *jj, ij, *esignal ;
double  t ;
{
    int     i, j, psignal ;
    double  diff, tolerance, meas_ang, cata_ang ;  
    void    mag_test() ;

    tolerance = TOL6 ; /* the tol. of difference in
                                 interangles between cat. and mesur. */

    meas_ang =  acos( m_star[n]->L[0] * m_star[m]->L[0]
	     +  m_star[n]->L[1] * m_star[m]->L[1]
	     +  m_star[n]->L[2] * m_star[m]->L[2])*180.0/PI*3600.0 ;
                        /*  the cosine of interangle  */
    
    for(i=0;i<ij;i++)
       for(j=i+1;j<ij;j++)
          {
          cata_ang =  acos( BLI[i].L[0] * BLI[j].L[0]
                    + BLI[i].L[1] * BLI[j].L[1]
                    + BLI[i].L[2] * BLI[j].L[2])*180.0/PI*3600.0 ;
                            /* of catalog star pairs */
	  diff = fabs(meas_ang - cata_ang) ;
                         /* compare the interangles */

          if(diff < tolerance*180.0/PI*3600.0) {
		mag_test(n,m,i,j,&psignal) ;
		if (psignal <= 2) 
                   { 
                   *ii = i ; *jj = j ;
                   *esignal = psignal  ; 
                   return ; 
                   }   /* mag_test was successful(2) or ambiguous(0) */
             } 
          } 
       *ii = i ;
       *jj = j ;
       *esignal = 3 ;  /* no match found */
       return ;
}

void compare1(t, n, m, ii, jj, ij, esignal)
int  n, m, ij, *ii, *jj, *esignal ;
double t  ;
{
    int     i, j, ji, psignal ;
    double  diff, tolerance, meas_ang, cata_ang ; 
    void    mag_test() ;

    tolerance = TOL6 ; /* the tol. of difference in
                          interangles between cat. and mesur. */
    meas_ang =  acos( m_star[n]->L[0] * m_star[m]->L[0]
	     +  m_star[n]->L[1] * m_star[m]->L[1]
	     +  m_star[n]->L[2] * m_star[m]->L[2])*180.0/PI*3600.0 ;
                        /*  the cosine of interangle  */
    
    for(i=*ii;i<ij;i++) {
       if(i > *ii)  ji = i+1   ;
       else         ji = *jj+1 ;        

    for(j=ji;j<ij;j++)
       {
       cata_ang =  acos( BLI[i].L[0] * BLI[j].L[0]
                + BLI[i].L[1] * BLI[j].L[1]
                + BLI[i].L[2] * BLI[j].L[2])*180.0/PI*3600.0 ;
                     /* of catalog star pairs */
       diff = fabs(meas_ang - cata_ang) ;
                      /* compare the interangles */

       if(diff < tolerance*180.0/PI*3600.0 ) {
             mag_test(n,m,i,j,&psignal) ;
	     if(psignal <= 2) { *ii = i ; *jj = j ; 
                                *esignal = psignal ;
                                return ; }
               } 
             }
       } 
       *ii = i ;   /* ij */
       *jj = j ;   /* ij */
       *esignal = 3 ;
       return ;
}

void compare2(n, m, jj, ij, sign, tolerance, nstar)
int     n, m, *jj, ij, *sign ;
double  tolerance ;
{
    int     j ;
    double  diff, mag_diff, meas_ang, cata_ang ; 

    meas_ang =  acos( m_star[n]->L[0] * m_star[m]->L[0]
	     +  m_star[n]->L[1] * m_star[m]->L[1]
	     +  m_star[n]->L[2] * m_star[m]->L[2])*180.0/PI*3600.0 ;
                        /*  the cosine of interangle  */
    for(j=0;j<ij;j++) {
      if(NEWBLI[0].id != j && NEWBLI[1].id != j 
          && NEWBLI[2].id != j && NEWBLI[3].id != j 
             && NEWBLI[4].id != j && NEWBLI[5].id != j)  
        {
        cata_ang =  acos( BLI[NEWBLI[0].id].L[0] * BLI[j].L[0]
                 + BLI[NEWBLI[0].id].L[1] * BLI[j].L[1]
                 + BLI[NEWBLI[0].id].L[2] * BLI[j].L[2])*180.0/PI*3600.0 ;
			/* of catalog star pairs */
	diff = fabs(meas_ang - cata_ang) ; /* compare the interangles */

        if(diff < tolerance*180.0/PI*3600.0)   {
           mag_diff = fabs(BLI[j].mag - m_star[m]->mag) ;
           if(mag_diff < Limit)  {
                *jj = j ;
                *sign = 1 ;
	        return ;
	        }    
	     }   
           }    
        }
      *jj = j ;
      *sign = 2 ;  
      return ;
}

void compare3(n, m, jj, sign, tolerance)
int     n, m, jj, *sign ;
double  tolerance ;
{
   int     k ;
   double  diff, meas_ang, cata_ang ; 
    
   meas_ang =  acos( m_star[n]->L[0] * m_star[m]->L[0]
            +  m_star[n]->L[1] * m_star[m]->L[1]
            +  m_star[n]->L[2] * m_star[m]->L[2])*180.0/PI*3600.0 ;
                        /*  the cosine of interangle  */
   cata_ang =  acos( BLI[NEWBLI[1].id].L[0] * BLI[jj].L[0]
	    +  BLI[NEWBLI[1].id].L[1] * BLI[jj].L[1]
	    +  BLI[NEWBLI[1].id].L[2] * BLI[jj].L[2])*180.0/PI*3600.0 ;
		  /* of catalog star pairs */
   diff = fabs(meas_ang - cata_ang) ;
		  /* compare the interangles */
   if(diff < tolerance*180.0/PI*3600.0 ) {
	   for(k=0;k<3;k++)
		    NEWBLI[m].L[k] = BLI[jj].L[k] ;
	   NEWBLI[m].num = BLI[jj].starnum ;
	   NEWBLI[m].mag  = BLI[jj].mag ;
	   NEWBLI[m].id = jj ;
           return ;  
           }  

   *sign = 2 ;  /*  not confirmed */
   return ;
}

void compare4(n, m, jj, ij, sign, tolerance)
int     n, m, *jj, ij, *sign ;
double  tolerance ;
{
    int     j ;
    double  diff, mag_diff, meas_ang, cata_ang ; 

    meas_ang =  acos( m_star[n]->L[0] * m_star[m]->L[0]
	     +  m_star[n]->L[1] * m_star[m]->L[1]
	     +  m_star[n]->L[2] * m_star[m]->L[2])*180.0/PI*3600.0 ;
    for(j=(*jj)+1;j<ij;j++)  
      if(NEWBLI[0].id != j && NEWBLI[1].id != j 
          && NEWBLI[2].id != j && NEWBLI[3].id != j 
             && NEWBLI[4].id != j && NEWBLI[5].id != j)  
        {
        cata_ang =  acos( BLI[NEWBLI[0].id].L[0] * BLI[j].L[0]
                 +  BLI[NEWBLI[0].id].L[1] * BLI[j].L[1]
                 +  BLI[NEWBLI[0].id].L[2] * BLI[j].L[2])*180.0/PI*3600.0 ;
        diff = fabs(meas_ang - cata_ang) ;
			/* compare the interangles */
        if(diff < tolerance*180.0/PI*3600.0)  {
             mag_diff = fabs(BLI[j].mag - m_star[m]->mag) ;
             if(mag_diff < Limit)  {
                   *jj = j ;
                   *sign = 1 ;
	           return ;
                   }   
	     }    
         }    
     
     *jj = j ;
     *sign = 2 ;  /* not found */
     return ;
}

void  mag_test(n, m, i, j, psignal)
int  n, m, i, j, *psignal ;
{
    int     k ;
    double  meas0, meas1, cata0, cata1 ;

    meas0 = m_star[n]->mag ;
    meas1 = m_star[m]->mag ;
    cata0 = BLI[i].mag ;
    cata1 = BLI[j].mag ;

    if ( fabs(cata0-meas0) < Limit && fabs(cata1-meas1) < Limit) 
       {
       if ( fabs(cata0-meas1) < Limit && fabs(cata1-meas0) < Limit) 
          *psignal=0 ; 
       else
          *psignal=2 ;
       for(k=0;k<3;k++) {
          NEWBLI[n].L[k] = BLI[i].L[k] ;
          NEWBLI[m].L[k] = BLI[j].L[k] ;  }
       NEWBLI[n].num = BLI[i].starnum ;
       NEWBLI[m].num = BLI[j].starnum ;
       NEWBLI[n].mag   = BLI[i].mag ;
       NEWBLI[m].mag   = BLI[j].mag ;
       NEWBLI[n].id = i ;
       NEWBLI[m].id = j ;
       return ; 
       }

    if ( fabs(cata0-meas1) < Limit && fabs(cata1-meas0) < Limit) 
       {
       if ( fabs(cata0-meas0) < Limit && fabs(cata1-meas1) < Limit) 
          *psignal=0 ;  /* practically, this case wouldn't happen */
       else
          *psignal=2 ;
       for(k=0;k<3;k++) {
          NEWBLI[n].L[k] = BLI[j].L[k] ;
          NEWBLI[m].L[k] = BLI[i].L[k] ;  }
       NEWBLI[n].num = BLI[j].starnum ;
       NEWBLI[m].num = BLI[i].starnum ;
       NEWBLI[n].mag   = BLI[j].mag ;
       NEWBLI[m].mag   = BLI[i].mag ;
       NEWBLI[n].id = j ;
       NEWBLI[m].id = i ;
       return ; 
       }
               
    *psignal = 3 ;   /* if magnitude test is failed */
    return ; 
}

void cartesian(alpha, delta, L) /* Compute unit vector of star position */
double  alpha, delta, L[3] ; 
{
    L[0] = cos(delta) * cos(alpha) ; /* usually on x axis */
    L[1] = cos(delta) * sin(alpha) ; /* usually on y axis */ 
    L[2] = sin(delta) ;  /* usually on z axis */ 
}

void obs_ang(nstar, x, y, xmag) 
int nstar ;
double x[NUMs], y[NUMs], xmag[NUMs] ;
{
   int     i, j ;
   double  x0[3], mag0, xy ;

   for(i=0;i<3;i++)  {
      x0[i] = m_star[0]->L[i] ; 
      for(j=1;j<nstar;j++) 
         m_star[j-1]->L[i] = m_star[j]->L[i] ; 
      m_star[nstar-1]->L[i] = x0[i] ; 
      } 
   mag0 = m_star[0]->mag ;
   for(j=1;j<nstar;j++)
      m_star[j-1]->mag = m_star[j]->mag ;
   m_star[nstar-1]->mag = mag0 ;

   xy = x[0] ; 
   for(j=1;j<nstar;j++) 
      x[j-1] = x[j] ; 
   x[nstar-1] = xy ; 

   xy = y[0] ; 
   for(j=1;j<nstar;j++) 
      y[j-1] = y[j] ; 
   y[nstar-1] = xy ; 

   xy = xmag[0] ; 
   for(j=1;j<nstar;j++) 
      xmag[j-1] = xmag[j] ; 
   xmag[nstar-1] = xy ; 
}  

void rearrange(nstar, x, y, xmag) 
int nstar ;
double x[NUMs], y[NUMs], xmag[NUMs] ;
{
   int     i ;
   double  x0[3], mag0, xy ;

   for(i=0;i<3;i++)  {
      x0[i] = m_star[2]->L[i] ; 
      m_star[2]->L[i] = m_star[5]->L[i] ; 
      m_star[5]->L[i] = x0[i] ; 
      } 
   mag0 = m_star[2]->mag ;
   m_star[2]->mag = m_star[5]->mag ;
   m_star[5]->mag = mag0 ;

   xy = x[2] ; 
   x[2] = x[5] ; 
   x[5] = xy ; 

   xy = y[2] ; 
   y[2] = y[5] ; 
   y[5] = xy ; 

   xy = xmag[2] ; 
   xmag[2] = xmag[5] ; 
   xmag[5] = xy ; 
}

int sym_test(dist_tol)
double dist_tol ;
{
   int     L=0, M=1, N=2 ;
   double  cata_ang[3] ;

   cata_ang[0] =  NEWBLI[L].L[0] * NEWBLI[M].L[0]
               +  NEWBLI[L].L[1] * NEWBLI[M].L[1]
               +  NEWBLI[L].L[2] * NEWBLI[M].L[2] ;
   if(cata_ang[0] >  1.0) cata_ang[0] =  1.0 ;
   if(cata_ang[0] < -1.0) cata_ang[0] = -1.0 ; 
   cata_ang[0] =  acos(cata_ang[0]) ;  
   cata_ang[1] =  NEWBLI[L].L[0] * NEWBLI[N].L[0]
               +  NEWBLI[L].L[1] * NEWBLI[N].L[1]
               +  NEWBLI[L].L[2] * NEWBLI[N].L[2] ;
   if(cata_ang[1] >  1.0) cata_ang[1] =  1.0 ;
   if(cata_ang[1] < -1.0) cata_ang[1] = -1.0 ; 
   cata_ang[1] =  acos(cata_ang[1]) ;  
   cata_ang[2] =  NEWBLI[M].L[0] * NEWBLI[N].L[0]
               +  NEWBLI[M].L[1] * NEWBLI[N].L[1]
               +  NEWBLI[M].L[2] * NEWBLI[N].L[2] ;
   if(cata_ang[2] >  1.0) cata_ang[2] =  1.0 ;
   if(cata_ang[2] < -1.0) cata_ang[2] = -1.0 ; 
   cata_ang[2] =  acos(cata_ang[2]) ;  

   if ( (fabs(cata_ang[0]-cata_ang[1]) < dist_tol) ||
        (fabs(cata_ang[0]-cata_ang[2]) < dist_tol) ||
        (fabs(cata_ang[1]-cata_ang[2]) < dist_tol) )
                   return(3) ; 

   return(0) ; 
} 
  
