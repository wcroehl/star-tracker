/* 11(TI), 12(TP), 13(TD), 21(G1), 22(G2), 23(G3) */
#include "./atti06.h"
#include "./miscell.c"
#include "./starID_pm.c"
#include "./starID_dm.c"
#include "./q_method.c"
#include "./cuvf06.c"
#define DEL_TIME  1000.0
#define FNLEN     256

main (int argc, char *argv[])
{
  FILE     *gyro, *sttr, *outmag, *outi, *outo, *outr, *outn,
           *outpos, *outsc, *rate, *outdist, *outgbias,
           *out1a, *out1b, *out2, *out3, *outM, *angvel,
           *gbias, *outds, *outdq, *inival, *qc, *merge, *ixy, *t_diff,
           *outPT1, *outPT2,*outPT3, *outes1, *outes2; /*add 6 files in this line for test, 7/28/23*/
  char     ID_loca, cday[5], cq[4][12], cb[3][8], cit[17], ciq[4][17], cgt[16],
           fn_dir[FNLEN], fn_in3[FNLEN], fn_in4[FNLEN], fn_in6[FNLEN],
           fn_tmp[FNLEN], fn_dir1[FNLEN], fn_dir2[FNLEN] ;
  int      i, j, k, nsig, rg, error, rtn,
           nstar, crf_count, *id_star, *id_body,
           read_CCD(), read_CCD0(), no_stars(), n_add, n_g, gread = 0,
	   read_gyro(), consec, sign, opp_sign, cnt,
           rmsnum = 0, rmslnum = 0, rmsznum = 0, rmszlnum = 0,
           QC1, QC2,
           cnt0, cuvfout = 1, ready_merge = 0,
           cuvf(), q_method(),
           c_dm = 0, c_dm_ided = 0, c_pm = 0, c_pm_ided = 0, id_ ;
  long     znum, numlines, b_est_cnt=0 ;
  double   ccd_time, ccd_time0, bore_body[3],
           cctime, t_old, t_elapse,
           pseudo_time,
           q[4], q0[4], qq[4], qq_old[4], qp[4], qp_in[4],
           M[3][3], A[3][3], P_theta[3][3], P[6][6],
           ccdstep, temp,
           vtime, sc_ra, sc_dec, RK,
           t_gyro, gt_old, g_angle[4], w[3], w0[3], u[3], b_add[3],
	   gb_accu[Numg][3], b_est[3], b_est_average[3], b_est_sum[3],
           zsum[3], z_bar[3], z[3],
           t_rms, rmssum[3], rmslsum[3], rmszsum[3], rmszscalar,
           rmszlsum[3], rmszlscalar, d_i,
           total_mstar, total_ided,
           OneSig1, OneSig2, OneSig3, atti_rms, dt_merge,
	   t0, t1, wallclock,
	   x[Ns], y[Ns], xmag[Ns] ;
  void     star(), celladj(), starcell(), starcell2() ;
  time_t   now, now2 ;              /* cpu time calculation */

  time(&now) ;
  fprintf(stderr, "It's now %s \n", ctime(&now)) ;
  t0 = clock() ;
  QC1=0;  /* add code 7/14/23 */

  strcpy (fn_dir, argv[1]);
  strcpy (fn_dir1, argv[2]);
  strcpy (fn_dir2, argv[3]);

  /* INPUT : GYRO */
  strcpy(fn_tmp, fn_dir1) ; strcat(fn_tmp, "/g_rate_lrs.dat") ;
  gyro = fopen(fn_tmp, "r");

  /* INPUT : LRS STAR */
  strcpy(fn_tmp, fn_dir1) ; strcat(fn_tmp, "/FOVs.dat") ;
  sttr = fopen(fn_tmp, "r");

  /* INPUT : INITIAL ATTITUDE */
  strcpy(fn_tmp, fn_dir1) ; strcat(fn_tmp, "/q_initial.txt") ;
  inival = fopen(fn_tmp, "r");

  /* OUTPUTS */
  strcpy(fn_tmp, fn_dir2) ; strcat(fn_tmp, "/inertials.dat") ;
  outi = fopen(fn_tmp, "w");

  strcpy(fn_tmp, fn_dir2) ; strcat(fn_tmp, "/rotatings.dat") ;
  outr = fopen(fn_tmp, "w");

  strcpy(fn_tmp, fn_dir2) ; strcat(fn_tmp, "/gyrorate.dat") ;
  rate = fopen(fn_tmp, "w");

  strcpy(fn_tmp, fn_dir2) ; strcat(fn_tmp, "/ckg_M.dat") ;
  outM = fopen(fn_tmp, "w");

  strcpy(fn_tmp, fn_dir2) ; strcat(fn_tmp, "/ckg_1a.dat") ;
  out1a = fopen(fn_tmp, "w");

  strcpy(fn_tmp, fn_dir2) ; strcat(fn_tmp, "/glas_pos.dat") ;
  outpos = fopen(fn_tmp, "w");

  strcpy(fn_tmp, fn_dir2) ; strcat(fn_tmp, "/ckg_1b.dat") ;
  out1b = fopen(fn_tmp, "w");

  strcpy(fn_tmp, fn_dir2) ; strcat(fn_tmp, "/observed.dat") ;
  outo = fopen(fn_tmp, "w");

  strcpy(fn_tmp, fn_dir2) ; strcat(fn_tmp, "/ckg_2.dat") ;
  out2 = fopen(fn_tmp, "w");

  strcpy(fn_tmp, fn_dir2) ; strcat(fn_tmp, "/rmsstar.dat") ;
  outds = fopen(fn_tmp, "w");

  strcpy(fn_tmp, fn_dir2) ; strcat(fn_tmp, "/rmsckg.dat") ;
  outdq = fopen(fn_tmp, "w");

  strcpy(fn_tmp, fn_dir2) ;    strcat(fn_tmp, "/IDed_num.dat") ;
  outn = fopen(fn_tmp, "w");

  strcpy(fn_tmp, fn_dir2) ;    strcat(fn_tmp, "/ckg_3.dat") ;
  out3 = fopen(fn_tmp, "w");

  strcpy(fn_tmp, fn_dir2) ;    strcat(fn_tmp, "/ckg_gbias.dat") ;
  gbias = fopen(fn_tmp, "w");

  strcpy(fn_tmp, fn_dir2) ;    strcat(fn_tmp, "/gbias.average") ;
  outgbias = fopen(fn_tmp, "w");

  strcpy(fn_tmp, fn_dir2) ;    strcat(fn_tmp, "/ckg_angvel.dat") ;
  angvel = fopen(fn_tmp, "w");

  strcpy(fn_tmp, fn_dir2) ;    strcat(fn_tmp, "/star_mag.dat") ;
  outmag = fopen(fn_tmp, "w");

  strcpy(fn_tmp, fn_dir2) ;    strcat(fn_tmp, "/dm_dist.dat") ;
  outdist = fopen(fn_tmp, "w");

  strcpy(fn_tmp, fn_dir2) ;    strcat(fn_tmp, "/ist_qc.dat") ;
  qc = fopen(fn_tmp, "w");

  strcpy(fn_tmp, fn_dir2) ;    strcat(fn_tmp, "/merge.dat") ;
  merge = fopen(fn_tmp, "w");

  strcpy(fn_tmp, fn_dir2) ;    strcat(fn_tmp, "/ided_xy.dat") ;
  ixy = fopen(fn_tmp, "w");

  strcpy(fn_tmp, fn_dir2) ;    strcat(fn_tmp, "/t_diff.dat") ;
  t_diff = fopen(fn_tmp, "w");


  /*added from line 141 to line 160, 7/28/23*/

  /* Test OUTPUTS */
  /*P_theta(P) Test OutPut*/
  strcpy(fn_tmp, fn_dir2); strcat(fn_tmp, "/outpt1.dat");
  outPT1 = fopen(fn_tmp, "w");

  strcpy(fn_tmp, fn_dir2); strcat(fn_tmp, "/outpt2.dat");
  outPT2= fopen(fn_tmp, "w");

  strcpy(fn_tmp, fn_dir2); strcat(fn_tmp, "/outpt3.dat");
  outPT3 = fopen(fn_tmp, "w");

  /*esignal Test Output*/
  strcpy(fn_tmp, fn_dir2); strcat(fn_tmp, "/outes1.dat");
  outes1 = fopen(fn_tmp, "w");

  strcpy(fn_tmp, fn_dir2); strcat(fn_tmp, "/outes2.dat");
  outes2 = fopen(fn_tmp, "w");

  /*****************************************************************/

  /* load all data in datafile in memory area */
  star(fn_dir) ;
  celladj(fn_dir) ;
  starcell(fn_dir) ;
  starcell2(fn_dir) ;

  id_star = (int *) malloc(Ns*sizeof(int)) ;
  if ( id_star == NULL) {
     fprintf(stdout, " FAIL : malloc 1 \n") ; exit(5) ; }
  id_body = (int *) malloc(Ns*sizeof(int)) ;
  if ( id_body == NULL) {
     fprintf(stdout, " FAIL : malloc 2 \n") ; exit(5) ; }
  m_star = (struct VecStar2 **) malloc(Ns*sizeof(struct VecStar2 *)) ;
  if ( m_star == NULL) {
     fprintf(stdout, " FAIL : malloc 3 \n") ; exit(5) ; }
  for(i=0;i<Ns;i++) {
     m_star[i] = (struct VecStar2 *) malloc(sizeof(struct VecStar2)) ;
     if ( m_star[i] == NULL) {
        fprintf(stdout, " FAIL : malloc (read_CCD) \n") ;
        exit(5) ; }}
  crf    = (struct VecStar **) malloc(Ns * sizeof(struct VecStar *)) ;
  if ( crf == NULL) {
     fprintf(stdout, " FAIL : malloc 4 \n") ; exit(5) ; }
  for(i=0;i<Ns;i++) {
     crf[i] = (struct VecStar *) malloc(sizeof(struct VecStar)) ;
     if ( crf[i] == NULL) {
        fprintf(stdout, " FAIL : malloc 4a \n") ; exit(5) ; }}
  c_star = (struct VecStar2 **) malloc(Ns * sizeof(struct VecStar2 *)) ;
  if ( c_star == NULL) {
     fprintf(stdout, " FAIL : malloc 5 \n") ; exit(5) ; }
  for(i=0;i<Ns;i++) {
     c_star[i] = (struct VecStar2 *) malloc(sizeof(struct VecStar2)) ;
     if ( c_star[i] == NULL) {
        fprintf(stdout, " FAIL : malloc 6 \n") ; exit(5) ; }}
  b_star = (struct VecStar2 **) malloc(Ns * sizeof(struct VecStar2 *)) ;
  if ( b_star == NULL) {
     fprintf(stdout, " FAIL : malloc 7 \n") ; exit(5) ; }
  for(i=0;i<Ns;i++) {
     b_star[i] = (struct VecStar2 *) malloc(sizeof(struct VecStar2)) ;
     if ( b_star[i] == NULL) {
        fprintf(stdout, " FAIL : malloc 8 \n") ; exit(5) ; }}

  id_star = (int *) malloc(Ns*sizeof(int)) ;
  for(i=0;i<3;i++) bore_body[i] = T_B[i][2] ;  /* (0,0,1) in the ST frame */

  ccdstep = T_STEP_CCD ;
  for(i=0;i<3;i++)  rmssum[i] = rmslsum[i] = 0.0 ;
  for(i=0;i<3;i++)  b_add[i] = 0.0 ;
  n_add = 0 ;
  n_g = 0 ;
  consec = 0 ;
  nsig = 0 ;
  rg = 0 ;
  crf_count = cnt = cnt0 = 0 ;
  ID_loca = 'F' ;  id_ = 99 ; /* Failed as initial: Used actually? */
  total_mstar = total_ided = 0.0 ;
  for(i=0;i<4;i++) g_angle[i] = 0.0 ;
  gt_old = 0.0 ;

  /* KF (t < KF) : in opt. bench coordinate frame */
  w0[2] = w0[1] = 1.e-4 ; w0[0] = 2.0 * PI / T /10.0 + 1.e-4;

  fscanf(inival, "%s", cit) ; /* not needed so far */
  for (i=0;i<4;i++) { fscanf(inival, "%s", ciq[i]) ;
		      q[i] = atof(ciq[i]) ; }

  b_est_average[0] = 0.01/3600.*PI/180. ;
  b_est_average[1] = 0.01/3600.*PI/180. ;
  b_est_average[2] = 0.01/3600.*PI/180. ;

  temp = sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3]) ;
  for(i=0;i<4;i++) q[i] /= temp ;
  for(i=0;i<4;i++) q0[i] = qq_old[i] = q[i] ;

  for(i=0;i<6;i++) for(j=0;j<6;j++) P[i][j] = 0.0 ;
  for(i=0;i<3;i++) P[i][i] = 1.0e-6 ;
  for(i=3;i<6;i++) P[i][i] = 1.0e-8 ;
  /* for(i=0;i<6;i++) P[i][i] = 1.0 ; */

  /* READ STAR DATA */
  ccd_time = -2000.0 ;
  do {
     nstar = read_CCD0(&ccd_time, sttr, x, y, xmag) ;
     } while (nstar < 3); /* || ccd_time < 2500.0) ; */

  /* READ GYRO DATA */
  for (i=0;i<1;i++) {
     gread = read_gyro(&t_gyro, u, b_est_average, w, gyro) ;
     if (gread == EOF) {
        fprintf(stdout, " No more gyro data\n") ;
        exit(3) ; }
     }

  while ( (t_gyro - ccd_time) > 0.05) {
     nstar = read_CCD(&ccd_time, sttr, x, y, xmag)  ;
     if (nstar > Ns)  { fprintf(stdout,"Too many stars at t=%8.2f\n", ccd_time) ;
                        exit(3) ; }
     else if (nstar == -999) {
        fprintf(stdout,"Empty IST data file before 'while' loop\n") ;
        exit(3) ; }
     else {
        total_mstar += (double) nstar ; }
     }

  while ( (ccd_time - t_gyro) > 0.05) {
     gread = read_gyro(&t_gyro, u, b_est_average, w, gyro) ;
     if (gread == 1) {
        fprintf(stdout, " Empty gyro data file at ccd_time = %8.2f\n", ccd_time) ;
        exit(3) ; }
     }

  ccd_time0 = ccd_time ;
  cctime = ccd_time + Ctime ; /* set up the pm range */
  t_rms = t_gyro ;
  if (ccd_time < KF)  for(i=0;i<3;i++) w[i] = w0[i];

  while( ccd_time <= TIMELIMIT )
    {
    d_i = 0.0 ;

    fprintf(t_diff, "%15.6f  %15.6f  %15.6f\n", ccd_time, t_gyro, ccd_time-t_gyro) ;
    temp = sqrt(q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2]+q0[3]*q0[3]) ;
    for(i=0;i<4;i++) q0[i] /= temp ;

    cnt = 0 ;

    if ( ccd_time < cctime && nstar >= 3)  {
       c_pm++ ;
       rtn = 0 ;
       if (nstar > 6) nstar = 6 ; /* up to 6 stars for PM */
       rtn = starID_pm(ccd_time, q0, nstar, &cnt, bore_body,
                 outi, outr, outo, outmag, ixy, x, y, xmag, outes1,outes2) ; /*added outes1 and outes2, 7/28/23*/

       /* fprintf(stderr, "P: ccd_time=%15.6f nstar=%3d  cnt = %3d\n", ccd_time, nstar, cnt) ; */

       for(i=0;i<Ns;i++) x[i] = y[i] = xmag[i] = 0.0 ;
       if (rtn != 1)
          {
          fprintf(stdout,"starID_pm() finished abnormally at %10.2f\n",
                                 ccd_time) ;
          exit(6) ;
          }

       if (cnt > 0) c_pm_ided++ ;
       crf_count = cnt ;

       ID_loca = 'P' ; id_ = 12 ;
       if (cnt >= 3 && rg > 0 )  /* UPDATE from > 0 */
                    {
                    cnt0 = n_add = n_g = rg = 0 ;
                    for(i=0;i<3;i++) b_add[i] = 0.0 ;
                    }
       }
    else if (ccd_time >= cctime && nstar > 0 && consec != 0)  {
       c_dm++ ;
       rtn = 0 ;
       rtn = starID_dm(ccd_time, q0, nstar, &cnt,
                     cctime, w, crf_count, id_star, id_body,
                     outi, outr, outo, outmag, outdist, ixy, x, y, xmag, 677, 700) ;

       for(i=0;i<Ns;i++) x[i] = y[i] = xmag[i] = 0.0 ;

       if (rtn != 1)
          {
          fprintf(stdout,"starID_dm() finished abnormally at %10.2f\n",
                                 ccd_time) ;
          exit(6) ;
          }

       if (cnt > 0) c_dm_ided++ ;
       crf_count = cnt ;
       ID_loca = 'D' ; id_ = 13 ;
       if (cnt > 1) {    /* UPDATE from 0 */
                    cnt0 = 0 ;
                    if (rg > 0)
                       {
                       n_add = n_g = rg = 0 ;
                       for(i=0;i<3;i++) b_add[i] = 0.0 ;
                       }
                    }
       }
    else { cnt = 0 ; } ;

    total_ided += (double) crf_count ;  /* just for statistics */

    if (cnt >= 3)    /* QUEST calculation for QC */
         {
        fprintf(outPT1, "%15.6f \n%20.15f %20.15f %20.15f \n%20.15f %20.15f %20.15f \n% 20.15f %20.15f %20.15f \n", ccd_time, P_theta[0][0], P_theta[0][1], P_theta[0][2], P_theta[1][0], P_theta[1][1], \
            P_theta[1][2], P_theta[2][0], P_theta[2][1], P_theta[2][2]);
        /*added the above fprintf function, 7/28/23*/

        rtn = q_method(ccd_time, qq, P_theta, cnt, &d_i, out1a, out2, outPT2, outPT3) ; /*added outPTs, 7/28/23*/
         if (rtn != 0)
            {
            fprintf(stdout,"q_method() finished abnormally at %10.2f\n",
                                         ccd_time) ;
            exit(6) ;
            }
         /* Discontinuity Check */
	 for (i=0;i<4;i++) if (qq[i] > 0.3 && qq[i]*qq_old[i] < 0) opp_sign++ ;
         if (opp_sign >= 2) for(i=0;i<4;i++)  sign = 1 ;
         else                                 sign = 0 ;

         opp_sign = 0 ;

         if (sign == 1)
            {
            if( nsig == 0)   nsig = 1 ;
            else             nsig = 0 ;
            sign = 0 ;
            }
         for(i=0;i<4;i++) qq_old[i] = qq[i] ;
         if( nsig == 1 ) for(i=0;i<4;i++)  qq[i] = -qq[i] ;

         fprintf(out1a,"%15.6f", ccd_time) ;
	 for(i=0;i<4;i++)
	    fprintf(out1a,"  %15.10f",qq[i]) ;
	 fprintf(out1a,"\n") ;

	 fprintf(out2,"%15.6f", ccd_time) ;
	 fprintf(out2,"   %12.8f    %12.8f    %12.8f ", \
	        sqrt(4*P_theta[0][0])*180./PI*3600,  \
	        sqrt(4*P_theta[1][1])*180./PI*3600,  \
	        sqrt(4*P_theta[2][2])*180./PI*3600 ) ;
	 fprintf(out2,"\n") ;

         if (consec < 9)  {
            consec++ ;
            if (consec == 9) { /* 10th consecutive data: (re)start filter */
               /* for(i=0;i<4;i++) qp[i] = q[i] = qq[i] ;  */
               for(i=0;i<6;i++) for(j=0;j<6;j++) P[i][j] = 0.0 ;
               for(i=0;i<3;i++) for(j=0;j<3;j++) P[i][j] = P_theta[i][j] ;
               t_old = ccd_time ;
               }

            q_to_A(qq, A) ;
            for(i=0;i<4;i++) qp[i] = q[i] = q0[i] = qq[i] ;
            for(i=0;i<3;i++) for(j=0;j<3;j++) M[i][j] = 0.0 ;
            for(i=0;i<3;i++) for(j=0;j<3;j++) for(k=0;k<3;k++)
                   M[i][j] += T_B[k][i] * A[k][j]   ;

            fprintf(out1b,"%15.6f  11 ", ccd_time) ;
            for(i=0;i<4;i++)  fprintf(out1b,"  %15.12f",qq[i]) ;
            fprintf(out1b,"\n") ;
            fprintf(qc,"%15.6f  %3d %6.2f\n",ccd_time, 1, -0.5) ;

            fprintf(outn,"%15.6f  11  %2d  %2d  %2d  %8.3f\n",
            ccd_time, nstar, cnt, consec, 0.0) ;

            fprintf(outM,"%12.6f %3d\n%20.15f %20.15f %20.15f \n%20.15f %20.15f %20.15f \n% 20.15f %20.15f %20.15f \n",
                  ccd_time, QC1, M[0][0], M[0][1], M[0][2], M[1][0], M[1][1],  \
                      M[1][2], M[2][0], M[2][1], M[2][2]) ;
            }
         }

    if (consec == 10 && cnt != 0) {           /* EKF */
         t_elapse = ccd_time - ccd_time0 ;    /* Tuning adjustment */

         for(i=0;i<3;i++) b_est[i] = b_est_average[i] ;
         cuvfout = cuvf(ccd_time, &t_old, t_elapse, qq, q, P, w, b_est,
                        cnt, zsum, &znum, ccdstep) ;

         for(i=0;i<3;i++) b_est_average[i] = b_est[i] ;
	 if (ccd_time >= 0.0)
            {
	    for(i=0;i<3;i++) b_est_sum[i] += b_est[i] ;
            b_est_cnt++ ;
            }

         if (cuvfout != 0) {
            fprintf(stdout," cuvf() finished abnormally at %10.2f\n",
                                         ccd_time) ;
            exit (6) ; }

         if (ready_merge == 1) {
	         fprintf(merge," %15.6f  %15.6f\n",
            	    ccd_time-dt_merge/10.0, ccd_time-1.0) ;
	         ready_merge = 0 ;
		 }

         if (ccd_time < KF)  for(i=0;i<3;i++) w0[i] = w[i];

         for(i=0;i<4;i++)  qp[i] = q0[i] = q[i] ;
         q_to_A(qp,  A) ;

         for(i=0;i<3;i++) for(j=0;j<3;j++) M[i][j] = 0.0 ;
         for(i=0;i<3;i++) for(j=0;j<3;j++) for(k=0;k<3;k++)
                          M[i][j] += T_B[k][i] * A[k][j] ;

         if (ccd_time > cctime)  /* for QC */
            {
            for(i=0;i<3;i++) rmssum[i]  += zsum[i] ; rmsnum  += znum ;
            for(i=0;i<3;i++) rmslsum[i] += zsum[i] ; rmslnum += znum ;
            }

         if (cnt >= 3)
            {
            for(i=0;i<3;i++) qp_in[i] = -qp[i] ;
                             qp_in[3] = qp[3] ;

            qcomp(qq,qp_in,z_bar) ;
            temp = sqrt(z_bar[0]*z_bar[0] + z_bar[1]*z_bar[1]
                 +      z_bar[2]*z_bar[2] + z_bar[3]*z_bar[3]) ;

            for(i=0;i<4;i++) z_bar[i] /= temp ;

            for(i=0;i<3;i++) z[i] = 2 * z_bar[i] ; /* 3 x 1 vector */

            if (t_gyro > cctime)   /* for QC */
               {
               for(i=0;i<3;i++) rmszsum[i] += sqrt(z[i]*z[i]) ;
               rmszscalar += sqrt(z_bar[3]*z_bar[3]) ;
               rmsznum++ ;
               for(i=0;i<3;i++) rmszlsum[i] += sqrt(z[i]*z[i]) ;
               rmszlscalar += sqrt(z_bar[3]*z_bar[3]) ;
               rmszlnum++ ;
               }
            }

         if ( cnt == 2 ) /* || (ID_loca == 'D' && cnt >= 3) ) */
            d_i = b_star[0]->L[0]*b_star[1]->L[0] +b_star[0]->L[1]*b_star[1]->L[1]
                + b_star[0]->L[2]*b_star[1]->L[2] ;
         if ( cnt  < 2 )  d_i = 0.0 ;

         fprintf(gbias, "%15.6f  %15.8e  %15.8e  %15.8e\n",
                 ccd_time, b_est[0]*3600/PI*180.0,
                           b_est[1]*3600/PI*180.0,
                           b_est[2]*3600/PI*180.0) ;
         fprintf(angvel, "%15.6f  %15.5f  %15.5f  %15.5f\n", \
                 ccd_time, w[0]*180.0/PI*3600.,
                           w[1]*180.0/PI*3600.,
                           w[2]*180.0/PI*3600.) ;

         OneSig1 = sqrt(P[0][0])*180./PI*3600 ;
         OneSig2 = sqrt(P[1][1])*180./PI*3600 ;
         OneSig3 = sqrt(P[2][2])*180./PI*3600 ;

         fprintf(out3,"%15.6f",ccd_time) ;
         fprintf(out3,"   %12.8f   %12.8f   %12.8f    %12.8f  %12.8f  %12.8f\n",
                 OneSig1, OneSig2, OneSig3,
                 sqrt(P[3][3])*180./PI*3600,  \
                 sqrt(P[4][4])*180./PI*3600,  \
                 sqrt(P[5][5])*180./PI*3600 ) ;

         atti_rms = sqrt(OneSig1*OneSig1+OneSig2*OneSig2) ;
         if   (atti_rms <= 2.0 )  QC1 = 2 ;
         else if (atti_rms > 2.0 && atti_rms < 4.0)  QC1 = 1 ;
         else QC1 = 0 ;

         if (ID_loca == 'P')  fprintf(out1b,"%15.6f  12 ",ccd_time) ;
         if (ID_loca == 'D')  fprintf(out1b,"%15.6f  13 ",ccd_time) ;
         for(i=0;i<4;i++)  fprintf(out1b,"  %15.12f",qp[i]) ;
         fprintf(out1b,"\n") ;
         fprintf(qc,"%15.6f  %3d %6.2f\n",ccd_time, QC1, atti_rms) ;

         fprintf(outM,"%12.6f %3d\n%20.15f %20.15f %20.15f \n%20.15f %20.15f %20.15f\n% 20.15f %20.15f %20.15f \n",
               ccd_time, QC1, M[0][0], M[0][1], M[0][2], M[1][0], M[1][1],  \
               M[1][2], M[2][0], M[2][1], M[2][2]) ;

         if (cnt > 2)
            fprintf(outn,"%15.6f  %2d  %2d  %2d  %2d  %8.3f\n",
               ccd_time, id_, nstar, crf_count,
               consec, acos(d_i/cnt)*180.0/PI);
         else if (cnt == 2)
            fprintf(outn,"%15.6f  %2d  %2d  %2d  %2d  %8.3f\n",
               ccd_time, id_, nstar, crf_count,
               consec, acos(d_i)*180.0/PI);
         else /* if (cnt == 1)  */
            fprintf(outn,"%15.6f  %2d  %2d  %2d  %2d  %8.3f\n",
               ccd_time, id_, nstar, crf_count,
               consec, 0.0);
         }  /* end of if ( consec == 10 && cnt != 0) */

    if (cnt == 0)
      {
      cnt0++ ;
      if (cnt0 > 100.0 || (cnt0 > 10.0 && nstar >= 3) ) {    /* Changed */
	  cctime = ccd_time + Ctime ;
	  consec = 0 ;
	  QC1 = 0 ; atti_rms = -1.0 ; }
      if (cnt0 > 600.0) {
	  /* fprintf(stdout, " 1. Warning: No obs for %6.2f min at %15.6f\n",
		cnt0/10.0/60.0, t_gyro) ; */
	  dt_merge = cnt0 ;
	  ready_merge = 1 ;
	  QC1 = 0 ; atti_rms = -2.0 ; }

      rg = no_stars(ccd_time, t_gyro, q0, w) ;

      for(i=0;i<4;i++) qp[i] = q0[i] ;
      q_to_A(qp, A) ;
      for(i=0;i<3;i++) for(j=0;j<3;j++) M[i][j] = 0.0 ;
      for(i=0;i<3;i++) for(j=0;j<3;j++) for(k=0;k<3;k++)
               M[i][j] += T_B[k][i] * A[k][j]   ;

      fprintf(out1b,"%15.6f  21 ", ccd_time) ; /* no star IDed */
      for(i=0;i<4;i++)  fprintf(out1b,"  %15.12f",qp[i]) ;
      fprintf(out1b,"\n") ;
      fprintf(qc,"%15.6f  %3d %6.2f\n",ccd_time, QC1, atti_rms) ;

      ID_loca = 'G' ; id_ = 21 ; consec = 0 ;
      fprintf(outn,"%15.6f  %2d  %2d  %2d  %2d  %8.3f\n",
             ccd_time, id_, nstar, cnt, consec, 0.0) ;

      fprintf(outM,"%12.6f %3d\n%20.15f %20.15f %20.15f \n%20.15f %20.15f %20.15f \n% 20.15f %20.15f %20.15f \n",
          ccd_time, QC1, M[0][0], M[0][1], M[0][2], M[1][0], M[1][1],  \
                  M[1][2], M[2][0], M[2][1], M[2][2]) ;
      }

    if (consec == 9 && t_old == ccd_time) consec = 10 ;

    while ( (nstar=read_CCD(&ccd_time, sttr, x, y, xmag)) == 0)
       {
       if ( nstar == -999 ) {
           fprintf(stdout,"No more IST data at or after %10.2f\n",ccd_time) ;
           if ( t_elapse < DEL_TIME )    exit(3) ;
           else                          break   ;
           }
       gread = read_gyro(&t_gyro, u, b_est_average, w, gyro) ;
       if (gread == 1) {
           fprintf(stdout,"No more gyro data. IST data may exist at %10.2f\n", t_gyro) ;
           if (t_gyro < 86400.0)  exit(3) ;
           else
   	     {
	     fprintf(outgbias, "%15ld  %15.5f  %15.5f  %15.5f\n",
			 b_est_cnt,  b_est_sum[0]/b_est_cnt*3600/PI*180.0,
	                 b_est_sum[1]/b_est_cnt*3600/PI*180.0,
		         b_est_sum[2]/b_est_cnt*3600/PI*180.0) ;
             exit(0) ;
	     }
           }
       if (ccd_time < KF)  for(i=0;i<3;i++) w[i] = w0[i];

       rg = no_stars(ccd_time, t_gyro, q0, w) ;

       for(i=0;i<4;i++) qp[i] = q0[i] ;

       q_to_A(qp, A) ;

       for(i=0;i<3;i++) for(j=0;j<3;j++) M[i][j] = 0.0 ;

       for(i=0;i<3;i++) for(j=0;j<3;j++) for(k=0;k<3;k++)
                M[i][j] += T_B[k][i] * A[k][j]   ;

       cnt0++ ;
       if (cnt0 > 100.0 || (cnt0 > 10.0 && nstar >= 3) ) {  /* Changed */
	        cctime = ccd_time + Ctime ;
                consec = 0 ;
	        QC1 = 0 ; atti_rms = -1.0 ; }
       if (cnt0 > 600.0) {
	        dt_merge = cnt0 ;
	        ready_merge = 1 ;
	        QC1 = 0 ; atti_rms = -2.0 ; }

       fprintf(out1b,"%15.6f  22 ",ccd_time) ; /* No star observed */
       for(i=0;i<4;i++)  fprintf(out1b,"  %15.12f",qp[i]) ;
       fprintf(out1b,"\n") ;

       ID_loca = '2' ;  id_ = 22 ; consec = 0 ;

       fprintf(qc,"%15.6f  %3d %6.2f\n",ccd_time, QC1, atti_rms) ;

       fprintf(outn,"%15.6f  %2d  %2d  %2d  %2d  %8.3f\n",
             ccd_time, id_, 0, 0, consec, 0.0) ;

       fprintf(outM,"%12.6f %3d\n%20.15f %20.15f %20.15f \n%20.15f %20.15f %20.15f \n% 20.15f %20.15f %20.15f \n",
        ccd_time, QC1, M[0][0], M[0][1], M[0][2], M[1][0], M[1][1], M[1][2],
                M[2][0], M[2][1], M[2][2]) ;

       } /* end of while (nstar = read_CCD(  ) ) */

    if (nstar > 0) total_mstar += (double) nstar ;

    if (nstar > Ns ) { fprintf(stdout,"Too many stars at t=%8.2f\n", ccd_time) ;
                       exit(3) ; }
    if (nstar == -999) break ; /* 2nd break : 1st break occurred with no CCD data */
    gread = read_gyro(&t_gyro, u, b_est_average, w, gyro) ;
    if (gread == 1) {
       fprintf(stdout, "No more gyro data at %10.2f while IST data remains\n", ccd_time) ;
       if (ccd_time >= TIMELIMIT)
	  {
          fprintf(outgbias, "%15ld  %15.5f  %15.5f  %15.5f\n",
	         b_est_cnt,  b_est_sum[0]/b_est_cnt*3600/PI*180.0,
	                     b_est_sum[1]/b_est_cnt*3600/PI*180.0,
	                     b_est_sum[2]/b_est_cnt*3600/PI*180.0) ;
	  exit(0) ;
	  }
       else                    exit(3) ;
       }

    while ( (ccd_time - t_gyro) > 0.0001)
       {

       pseudo_time = t_gyro ;

       if (pseudo_time < KF)  for(i=0;i<3;i++) w[i] = w0[i];
       rg = no_stars(pseudo_time, t_gyro, q0, w) ;

       for(i=0;i<4;i++) qp[i] = q0[i] ;

       q_to_A(qp, A) ;

       for(i=0;i<3;i++) for(j=0;j<3;j++) M[i][j] = 0.0 ;

       for(i=0;i<3;i++) for(j=0;j<3;j++) for(k=0;k<3;k++)
                M[i][j] += T_B[k][i] * A[k][j]   ;

       cnt0++ ;
       if (cnt0 > 100.0) {
	        cctime = pseudo_time + Ctime ;
                consec = 0 ;
	        QC1 = 0 ; atti_rms = -1.0 ; }
       if (cnt0 > 600.0) {
	        dt_merge = cnt0 ;
	        ready_merge = 1 ;
	        QC1 = 0 ; atti_rms = -2.0 ; }

       fprintf(out1b,"%15.6f  23 ", pseudo_time) ;  /* No IST timetags */
       for(i=0;i<4;i++)  fprintf(out1b,"  %15.12f",qp[i]) ;
       fprintf(out1b,"\n") ;

       ID_loca = '3' ; id_ = 23 ;

       fprintf(qc,"%15.6f  %3d %6.2f\n", pseudo_time, QC1, atti_rms) ;

       fprintf(outn,"%15.6f  %2d  %2d  %2d  %2d  %8.3f\n",
               pseudo_time, id_, 0, 0, consec, 0.0) ;

       fprintf(outM,"%12.6f %3d\n%20.15f %20.15f %20.15f \n%20.15f %20.15f %20.15f \n% 20.15f %20.15f %20.15f \n",
         pseudo_time, QC1, M[0][0], M[0][1], M[0][2], M[1][0], M[1][1], M[1][2],
                M[2][0], M[2][1], M[2][2]) ;

       gread = read_gyro(&t_gyro, u, b_est_average, w, gyro) ;
       if (gread == 1) {
           fprintf(stdout,"No more gyro data. IST data may exist at %10.2f\n", t_gyro) ;
	   if (ccd_time >= 86400.)
	      {
              fprintf(outgbias, "%15ld  %15.5f  %15.5f  %15.5f\n",
	            b_est_cnt,  b_est_sum[0]/b_est_cnt*3600/PI*180.0,
	                        b_est_sum[1]/b_est_cnt*3600/PI*180.0,
	                        b_est_sum[2]/b_est_cnt*3600/PI*180.0) ;
	      exit(0) ;
	      }
           else                    exit(3) ;
           }

       } /* end of while (ccd_time == t_gyro) */

    /* if (t_gyro < KF)  for(i=0;i<3;i++) w[i] = w0[i]; */

    if ( (t_gyro - t_rms) >= DT_RMS) /* rms calc. with 10 minute intervals */
       {
       if (rmslnum > 0)
          fprintf(outds," %12.5f %9d %12.8f %12.8f %12.8f\n",
            ccd_time, rmslnum,
            rmslsum[0]/rmslnum*180.0/PI*3600.0,
            rmslsum[1]/rmslnum*180.0/PI*3600.0,
            rmslsum[2]/rmslnum*180.0/PI*3600.0) ;
       else
          fprintf(outds," %12.5f %9d %12.8f %12.8f %12.8f\n",
            ccd_time, rmslnum, 0.0, 0.0, 0.0) ;
       for(i=0;i<3;i++) rmslsum[i] = 0.0 ;
       rmslnum = 0 ;
       if (rmszlnum > 0)
          fprintf(outdq," %12.5f %9d %12.8f %12.8f %12.8f %12.8f\n",
            ccd_time, rmszlnum,
            rmszlsum[0]/rmszlnum*180.0/PI*3600.0,
            rmszlsum[1]/rmszlnum*180.0/PI*3600.0,
            rmszlsum[2]/rmszlnum*180.0/PI*3600.0,
            rmszscalar/rmsznum ) ;
       else
          fprintf(outdq," %12.5f %9d %12.8f %12.8f %12.8f %12.8f\n",
            ccd_time, rmszlnum, 0.0, 0.0, 0.0, 0.0) ;
       for(i=0;i<3;i++) rmszlsum[i] = 0.0 ;
       rmszlscalar = 0.0 ;
       rmszlnum = 0 ;
       t_rms += DT_RMS ;
       }
    if(ccd_time < KF) for(i=0;i<3;i++) w0[i] = w[i] ;
    numlines++ ;
    }  /* end of while */

  if (rmslnum > 0)
     fprintf(outds," %12.5f %9d %12.8f %12.8f %12.8f\n",
        ccd_time, rmslnum,
        rmslsum[0]/rmslnum*180.0/PI*3600.0,
        rmslsum[1]/rmslnum*180.0/PI*3600.0,
        rmslsum[2]/rmslnum*180.0/PI*3600.0) ;
  else
     fprintf(outds," %12.5f %9d %12.8f %12.8f %12.8f\n",
        ccd_time, rmslnum, 0.0, 0.0, 0.0) ;
  fprintf(outds," %12.5f %9d %12.8f %12.8f %12.8f\n",
        ccd_time, rmsnum,
        rmssum[0]/rmsnum*180.0/PI*3600.0,
        rmssum[1]/rmsnum*180.0/PI*3600.0,
        rmssum[2]/rmsnum*180.0/PI*3600.0) ;
  if (rmszlnum > 0)
     fprintf(outdq," %12.5f %9d %12.8f %12.8f %12.8f %12.8f\n",
        ccd_time, rmszlnum,
        rmszlsum[0]/rmszlnum*180.0/PI*3600.0,
        rmszlsum[1]/rmszlnum*180.0/PI*3600.0,
        rmszlsum[2]/rmszlnum*180.0/PI*3600.0,
        rmszscalar/rmsznum ) ;
  else
     fprintf(outdq," %12.5f %9d %12.8f %12.8f %12.8f %12.8f\n",
        ccd_time, rmszlnum, 0.0, 0.0, 0.0, 0.0) ;
  fprintf(outdq," %12.5f %9d %12.8f %12.8f %12.8f %12.8f\n",
        ccd_time, rmsznum,
        rmszsum[0]/rmsznum*180.0/PI*3600.0,
        rmszsum[1]/rmsznum*180.0/PI*3600.0,
        rmszsum[2]/rmsznum*180.0/PI*3600.0,
        rmszscalar/rmsznum ) ;

  fprintf(outgbias, "%15ld  %15.5f  %15.5f  %15.5f\n",
	   b_est_cnt,  b_est_sum[0]/b_est_cnt*3600/PI*180.0,
	               b_est_sum[1]/b_est_cnt*3600/PI*180.0,
	               b_est_sum[2]/b_est_cnt*3600/PI*180.0) ;

  fprintf(stdout, "End of PROGRAM\n") ;

  fprintf(stdout, "total_mstar = %15.0f  total_ided = %15.0f\n",
              total_mstar, total_ided) ;
  fprintf(stdout, "c_pm = %8d  c_pm_ided = %8d  c_dm = %8d  c_dm_ided = %8d\n",
              c_pm, c_pm_ided, c_dm, c_dm_ided) ;

  t1 = clock() ;
  time (&now2) ;
  wallclock = difftime (now2, now) ;
  fprintf(stderr, "It's now %s \n", ctime(&now2)) ;
  fprintf(stderr, "it took %.5f CPU secs.\n",(float)(t1-t0)/CLOCKS_PER_SEC) ;
  fprintf(stderr, "it took %.5f wall clock secs.\n", wallclock) ;

  exit (0) ;
 } /* End of Program */

void star(fn_dir)
char fn_dir[FNLEN] ;
  {
  FILE  *in ;
  static char  catnum[7], catmom[6], rig_asc[13], decli[13], magni[8] ;
  char         fn_tmp[FNLEN] ;
  int   i ;

  strcpy(fn_tmp, fn_dir) ; strcat(fn_tmp, "/lib/lrs_5p25.cat") ;
  if ((in = fopen(fn_tmp,"r")) == (FILE *) NULL)
     {
     fprintf(stdout, "Can't open 'lrs_5p25.cat'\n") ;
     exit (1) ; }

  for(i=0 ; i < NUMSTAR ; i++) {
     fscanf(in,"%s %s %s %s %s", catnum, catmom, rig_asc, decli, magni);
     stars[i].cat_num = atoi(catnum) ;
     stars[i].cat_mom = atoi(catmom) ;
     stars[i].ra = atof(rig_asc) ;
     stars[i].dec = atof(decli) ;
     stars[i].mag = atof(magni) ;
     }
  fclose(in) ;
  }

void starcell(fn_dir)
  {
  FILE  *in ;
  char  cellid[5], numstars[5], starnum[5], fn_tmp[FNLEN];
  int   i, j ;

  strcpy(fn_tmp, fn_dir) ; strcat(fn_tmp, "/lib/scel525.dat") ;
  if ((in = fopen(fn_tmp,"r")) == (FILE *) NULL)
     {
     fprintf(stdout, "Can't open 'scel.dat'\n") ;
     exit (1) ; }

  for(i=0 ; i < 225 ; i++)
     {
     fscanf(in, "%s %s",cellid, numstars) ;
     scell[i].cell_id = atoi(cellid) ;
     scell[i].num_stars = atoi(numstars) ;
     scell[i].star_num =
          (int *)  malloc(sizeof(int)*(scell[i].num_stars));

     if ( scell[i].star_num == NULL) {
         fprintf(stdout, " FAIL : malloc %3d (starcell) \n", i) ;
         exit(5) ; }

     for(j=0 ; j < scell[i].num_stars ; j++)
        {
        fscanf(in,"%s",starnum) ;
        scell[i].star_num[j] = atoi(starnum) ;
        }
     }
  fclose(in) ;
  }

void starcell2(fn_dir)
  {
  FILE  *in ;
  char  cellid[3], numstars[8], starnum[7], fn_tmp[FNLEN];
  int   i, j ;

  strcpy(fn_tmp, fn_dir) ; strcat(fn_tmp, "/lib/x525_tab.dat") ;
  if ((in = fopen(fn_tmp,"r")) == (FILE *) NULL)
     {
     fprintf(stdout, "Can't open 'x525_tab.dat'\n") ;
     exit (1) ; }

  for(i=0 ; i < 34 ; i++)
     {
     fscanf(in, "%s %s", cellid, numstars) ;
     scell2[i].cell_id = atoi(cellid) ;
     scell2[i].num_stars = atoi(numstars) ;
     scell2[i].star_num =
              (int *) malloc(sizeof(int)*(scell2[i].num_stars));

     if ( scell2[i].star_num == NULL) {
        fprintf(stdout, " FAIL : malloc %3d (starcell2) \n", i) ;
        exit(5) ; }

     for(j=0 ; j < scell2[i].num_stars ; j++)
        {
        fscanf(in,"%s", starnum) ;
        scell2[i].star_num[j] = atoi(starnum) ;
        }
     }
 fclose(in) ;
 }

void celladj(fn_dir)
  {
  FILE   *in ;
  char   numcell[3], adjcells[5], fn_tmp[FNLEN];
  int    i, j ;

  strcpy(fn_tmp, fn_dir) ; strcat(fn_tmp, "/lib/sadj.dat") ;
  if ((in = fopen(fn_tmp,"r")) == (FILE *) NULL)
     {
     fprintf(stdout, "Can't open 'sadj.dat'\n") ;
     exit (1) ; }

  for(i=0 ; i < 225  ; i++)
     {
     fscanf(in, "%s ", numcell) ;
     adjcell[i].num_adj = atoi(numcell) ;
     adjcell[i].cell_num =
             (int *) malloc(sizeof(int)*adjcell[i].num_adj)  ;

     if ( adjcell[i].cell_num == NULL) {
        fprintf(stdout, " FAIL : malloc %3d (celladj) \n", i) ;
        exit(5) ; }

     for(j=0; j < adjcell[i].num_adj ; j++)
        {
        fscanf(in,"%s", adjcells) ;
        adjcell[i].cell_num[j] = atoi(adjcells) ;
        }
     }
  fclose(in) ;
  }

int read_CCD(t, sttr, x, y, xmag)
FILE     *sttr ;
double   *t, x[Ns], y[Ns], xmag[Ns] ;
 {
    int     i, j, nstar, sm ;
    char    ct[17], cx[Ns][13], cy[Ns][13], cm[Ns][7] ;
    static int  dup_cnt = 0, rev_cnt = 0 ;
    double  phi[Ns], lambda[Ns], mag[Ns], s[2] ;
    void    sort_CCD() ;

    sm = 0 ;  /* 0:no problem, 1:bad data, 2:no data */
    if ( fscanf(sttr, "%s", ct) != EOF ) {
       *t = atof(ct) ;
       for(i=0;i<Ns;i++) { fscanf(sttr, "%s %s %s", cx[i], cy[i], cm[i]) ;
                          mag[i] = atof(cm[i]) ;
                         }

       nstar = 0 ;

       for(i=0;i<Ns;i++) {
          for(j=0;j<3;j++) { m_star[i]->L[j] = 0.0 ;
                             m_star[i]->mag = mag[i] ; }
          x[i] = y[i] = 0.0 ;
	  xmag[i] = mag[i] ;
 	  phi[nstar] = atof(cx[i]) ;
	  lambda[nstar] = atof(cy[i]) ;
	  mag[nstar] = mag[i] ;
	  if (mag[nstar] >= -2.0 && mag[nstar] < 7.5)  nstar++ ;
	  }
       for(i=0;i<nstar;i++) {
	  x[i] = phi[i] ;
	  y[i] = lambda[i] ;
	  xmag[i] = mag[i] ;
          phi[i] = phi[i]*PI/180.0 ; /* degree to radian */
          lambda[i] = lambda[i]*PI/180.0 ; /* degree to radian */

          s[0] = phi[i] ;     /* atan on tan(phi[i]) ; */
          s[1] = lambda[i] ;  /* atan on tan(lambda[i]) ; */
          m_star[i]->L[2] = 1./sqrt(1+s[0]*s[0]+s[1]*s[1]) ;
          m_star[i]->L[0] = s[0]*m_star[i]->L[2] ;
	  m_star[i]->L[1] = s[1]*m_star[i]->L[2] ;

          m_star[i]->mag = mag[i] ;
          }  /* end of for(i) */
    }  /* end of if(fread) */
    else {
       fprintf(stdout,"No more CCD measurements (read_CCD:0).\n") ;
       exit(100) ;
       return(-999) ;
       }
    sort_CCD(nstar, x, y, xmag) ;
    return(nstar) ;
 }

int read_CCD0(t, sttr, x, y, xmag)
FILE     *sttr ;
double   *t, x[Ns], y[Ns], xmag[Ns] ;
 {
    int     i, j, nstar, sm ;
    char    ct[17], cx[Ns][13], cy[Ns][13], cm[Ns][7] ;
    static int  dup_cnt = 0, rev_cnt = 0 ;
    double  phi[Ns], lambda[Ns], mag[Ns], s[2] ;
    void    sort_CCD() ;

    nstar = 0 ;
    sm = 0 ;  /* 0:no problem, 1:bad data, 2:no data */
    if (fscanf(sttr, "%s", ct) != EOF )
       {
       *t = atof(ct) ;
       for(i=0;i<Ns;i++) {
		 fscanf(sttr, "%s %s %s", cx[i], cy[i], cm[i]) ;
                 mag[i] = atof(cm[i]) ;
		 }

       for(i=0;i<Ns;i++) {
          for(j=0;j<3;j++) { m_star[i]->L[j] = 0.0 ;
                             m_star[i]->mag = mag[i] ; }
          x[i] = y[i] = 0.0 ;
	  xmag[i] = mag[i] ;
 	  phi[nstar] = atof(cx[i]) ;
	  lambda[nstar] = atof(cy[i]) ;
	  mag[nstar] = mag[i] ;
	  if (mag[nstar] >= -2.0 && mag[nstar] < 7.5)  nstar++ ;
	  }
       for(i=0;i<nstar;i++) {
	  x[i] = phi[i] ;
	  y[i] = lambda[i] ;
	  xmag[i] = mag[i] ;
          phi[i] = phi[i]*PI/180.0 ; /* degree to radian */
          lambda[i] = lambda[i]*PI/180.0 ; /* degree to radian */

          s[0] = phi[i] ;    /* tan(phi) in original code */
          s[1] = lambda[i] ; /* tan(lambda) in original code */
          m_star[i]->L[2] = 1./sqrt(1+s[0]*s[0]+s[1]*s[1]) ;
          m_star[i]->L[0] = s[0]*m_star[i]->L[2] ;
	  m_star[i]->L[1] = s[1]*m_star[i]->L[2] ;

          m_star[i]->mag = mag[i] ;
          }  /* end of for(i) */
    }  /* end of if(fread) */
    else {
       fprintf(stdout,"No more CCD measurements (read_CCD:0).\n") ;
       return(-999) ;
       }
    sort_CCD(nstar, x, y, xmag) ;
    return(nstar) ;
 }

void sort_CCD(nstar, x, y, xmag)
int  nstar ;
double x[Ns], y[Ns], xmag[Ns] ;
 {
    int  i, j, k ;
    double tmp_x, tmp_y, tmp_mag ;
    struct VecStar2 *temporary ;

    temporary = (struct VecStar2 *) malloc(sizeof(struct VecStar2)) ;
    if ( temporary == NULL) {
       fprintf(stdout,"FAIL : malloc (sort_CCD)\n"); exit(5) ;}
    for(i=0;i<nstar;i++) {
       for(j=i+1;j<nstar;j++)  /* bright star first */
           if(m_star[j]->mag < m_star[i]->mag)
              {
              temporary->mag     = m_star[i]->mag     ;
              for(k=0;k<3;k++) temporary->L[k] = m_star[i]->L[k] ;
	      tmp_x = x[i] ; tmp_y = y[i] ; tmp_mag = xmag[i] ;

              m_star[i]->mag     = m_star[j]->mag     ;
              for(k=0;k<3;k++) m_star[i]->L[k] = m_star[j]->L[k] ;
	      x[i] = x[j] ; y[i] = y[j] ; xmag[i] = xmag[j] ;

              m_star[j]->mag     = temporary->mag     ;
              for(k=0;k<3;k++) m_star[j]->L[k] = temporary->L[k] ;
	      x[j] = tmp_x ; y[j] = tmp_y ; xmag[j] = tmp_mag ;
              }
       }
}

int read_gyro(t_gyro, u, b, w, gyro)
FILE      *gyro ;
double    *t_gyro, u[3], b[3], w[3] ;
  {
  char cgt[16], cx[16], cy[16], cz[16] ;
  int  i ;

  if (fscanf(gyro, "%s", cgt) != EOF)
    {
    *t_gyro = atof(cgt) ;
    fscanf(gyro, "%s %s %s", cx, cy, cz) ;

    u[0] = atof(cx)/3600*PI/180. ;
    u[1] = atof(cy)/3600*PI/180. ;
    u[2] = atof(cz)/3600*PI/180. ;
    }

  else { fprintf(stdout, " No gyro measurement (read_gyro).\n") ;
         return(1) ; }

  for(i=0;i<3;i++) w[i] = u[i] - b[i] ;
  return (0) ;
}

int no_stars(t, t_gyro, q, w)
double  t, t_gyro, q[4], w[3] ;
{
   int     i, j ;
   double  w_mag, C, S, M[4][4], q_new[4], dq[4], temp ;

      temp = sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3]) ;
      for(i=0;i<4;i++) q[i] /= temp ;

      {

      dq[0] =( w[2]*q[1]-w[1]*q[2]+w[0]*q[3])/2.0/10.0 ;
      dq[1] =(-w[2]*q[0]+w[0]*q[2]+w[1]*q[3])/2.0/10.0 ;
      dq[2] =( w[1]*q[0]-w[0]*q[1]+w[2]*q[3])/2.0/10.0 ;
      dq[3] =(-w[0]*q[0]-w[1]*q[1]-w[2]*q[2])/2.0/10.0 ;

      for(i=0;i<4;i++) q[i] += dq[i] ;

      temp = sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3]) ;
      for(i=0;i<4;i++) q[i] /= temp ;
      }

      temp = sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3]) ;
      for(i=0;i<4;i++) q[i] /= temp ;
   return 2 ;
}

