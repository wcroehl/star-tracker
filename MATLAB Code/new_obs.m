function [rtn] = new_obs(t, t_start, t_end, count, candi_num, oldnum,...
         new, new_ra, new_dec, new_mag, id_star, id_count, vici, Mp, outdist, scell2, stars)
% FILE   *outdist ;
% struct VecStar *new ;
% int    *count, *candi_num, oldnum, *id_star, id_count ;
% double t, t_start, t_end, *new_ra, *new_dec, *new_mag, vici, Mp[3][3] ; 
%  {
%  static int  previ_N_zone, thisnum ;
%  int    i, j, *intptr, N_zone, astaragain, thisnum0, thisnum1, number,
%         zone() ;
%  double BD_vec[3], BD_ra, BD_dec, 
%         ra_diff_new, ra_diff_BD, i_ra, i_dec, i_mag, i_L[3],
%         vec_dot, arc_dist ;
 rtn = -1;
 previ_N_zone = 0;
 astaragain = 0 ;
 thisnum = 0;
 thisnum0 = 0;
 BD_limit = 6;
 FOV_limit = 12;
 for i = 0:1:2 
    BD_vec(i+1) = Mp(2,i+1) ;
 end %/* give star vector in the CRF frame */ 
 BD_ra = atan2(BD_vec(2), BD_vec(1)) * 180.0 / pi ;
               if ( BD_ra  <  0.0)
                    BD_ra = BD_ra + 360.0 ;
               end
 BD_dec = asin(BD_vec(3)) * 180.0 / pi ;

 N_zone = zone(BD_dec) ; %/* N_zone_ON */

 if (N_zone == -1)
   fprintf(stdout,"find()::zone() didn't finish normally.\n") ;
   return;
 end
 if (N_zone ~= previ_N_zone)
     cell = scell2(N_zone-1).num_stars;
   for j=1:1:cell
       if (scell2(N_zone-1).star_num(j) == oldnum.num)
           thisnum = j+1 ;
       end
       if (t > t_start && t < t_end) 
           fprintf(stdout, "Here? j=%4d\n", j) ;
           break ;
       end
   end
 end

 previ_N_zone = N_zone ;

 intptr = candi_num ;

 if (thisnum == 0) 
     thisnum = 1 ;
 end
 if (thisnum > scell2(N_zone-1).num_stars) 
     thisnum = scell2(N_zone-1).num_stars ;  
     if (t > t_start && t < t_end) 
         fprintf(stdout, "Here2? thisnum=%4d\n", thisnum) ;
     end
 end
     
 thisnum1 = thisnum0;
 thisnum0 = thisnum;

 if ( BD_dec >= 0.0) 
     BD_dec = BD_dec + 6.0; %/* 4.0 */
 end
 if ( BD_dec <  0.0)
     BD_dec = BD_dec - 6.0; %/* 4.0 */
 end

 i=0;
 while (true)
     
     i_ra  = stars(scell2(N_zone-1).star_num(thisnum-i)).ra;
     i_dec = stars(scell2(N_zone-1).star_num(thisnum-i)).dec;
     i_mag = stars(scell2(N_zone-1).star_num(thisnum-i)).mag;
  
     ra_diff_new = abs(new_ra - i_ra);
     ra_diff_BD  = abs(BD_ra  - i_ra) ;
     if (ra_diff_new > 180.0)  
         ra_diff_new = 360.0 - ra_diff_new ;
     end
     if (ra_diff_BD  > 180.0)  
         ra_diff_BD  = 360.0 - ra_diff_BD  ;
     end
     ra_diff_new = ra_diff_new * cos(new_dec*pi/180.0) * 3600.0 ;
     ra_diff_BD = ra_diff_BD * cos(  BD_dec*pi/180.0) * 3600.0 ;

     if ( abs(new_dec) > 80.0 || abs(BD_dec) > 80.0)   %/* 88 */
         ra_diff_BD  = 100.0 ;
         ra_diff_new = 100.0 ;
     end

     i_L(1) = cos(i_dec*pi/180.0)*cos(i_ra*pi/180.0) ;
     i_L(2) = cos(i_dec*pi/180.0)*sin(i_ra*pi/180.0) ;
     i_L(3) = sin(i_dec*pi/180.0) ;

     vec_dot = new.L(1)*i_L(1) + new.L(2)*i_L(2) + new.L(3)*i_L(3) ;
     if (vec_dot > 1.0)  
         vec_dot =  1.0 ;
     end
     if (vec_dot < -1.0) 
         vec_dot = -1.0 ;
     end
     arc_dist = acos(vec_dot)*180.0/pi*3600.0 ;

     %/* number = scell2[N_zone-1].star_num[thisnum+i] ; */

     if(abs(arc_dist) < vici && abs(i_mag - new_mag) < mTOL) %/* a possible match */
        if (t > t_start && t < t_end) 
            fprintf(stdout, "Here+-? thisnum=%4d  i=%3d\n", thisnum, i);
            thisnum1 = thisnum-i+1;
        end
        for j=1:1:id_count
           if(id_star(j) == stars(scell2(N_zone-1).star_num(thisnum-i)).cat_num) 
             astaragain = 1;
           end
        end
        if (astaragain ~= 1)  %/* no repeated star */
           intptr  = stars(scell2(N_zone-1).star_num(thisnum-i)).cat_num;
           fprintf(outdist, "2  %15.6f %5d %15.8f \n", t, intptr, arc_dist) ;
           count = count + 1;
           intptr = intptr + 1;
        end
      end
      i = i + 1;
      astaragain = 0 ;
      if( (thisnum - i) < 0 )  %/* break ;  */
            thisnum = thisnum + scell2(N_zone-1).num_stars ;  
            if (t > t_start && t < t_end) 
                fprintf(stdout, "Here3? thisnum=%4d\n", thisnum) ;
            end
      end
      if (t > t_start && t < t_end)
        fprintf(stdout, "+++ ra_diff_BD = %9.2f, ra_diff_new = %9.2f\n",...
               ra_diff_BD, ra_diff_new);
        fprintf(stdout, "  i = %3d  <<<<  i_con = %6.2f\n",... 
           i, scell2(N_zone-1).num_stars/2.0) ;
      end
      if(~(abs(ra_diff_BD) < BD_limit*3600+1800.) || ...
         ~((abs(ra_diff_new) < FOV_limit*3600+1800.) &&  i < floor(scell2(N_zone-1).num_stars/2.0))) break; end
 end
 kkk=0;
  while(((abs(ra_diff_BD) < BD_limit*3600+1800.) || (abs(ra_diff_new) < FOV_limit*3600+1800.))... 
         &&  i < floor(scell2(N_zone-1).num_stars/2.0))
     
     i_ra  = stars(scell2(N_zone-1).star_num(thisnum-i)-1).ra  ;
     i_dec = stars(scell2(N_zone-1).star_num(thisnum-i)-1).dec ;
     i_mag = stars(scell2(N_zone-1).star_num(thisnum-i)-1).mag ;

     ra_diff_new = abs(new_ra - i_ra) ;
     ra_diff_BD  = abs(BD_ra  - i_ra) ;
     if (ra_diff_new > 180.0)  
         ra_diff_new = 360.0 - ra_diff_new ;
     end
     if (ra_diff_BD  > 180.0)  
         ra_diff_BD  = 360.0 - ra_diff_BD  ;
     end
     ra_diff_new = ra_diff_new * cos(new_dec*pi/180.0) * 3600.0 ;
     ra_diff_BD  = ra_diff_BD * cos(BD_dec*pi/180.0) * 3600.0 ;

     if (abs(new_dec) > 80.0 || abs(BD_dec) > 80.0)   %/* 88 */
        ra_diff_BD  = 100.0 ;
        ra_diff_new = 100.0 ;
     end

     i_L(1) = cos(i_dec*pi/180.0)*cos(i_ra*pi/180.0) ;
     i_L(2) = cos(i_dec*pi/180.0)*sin(i_ra*pi/180.0) ;
     i_L(3) = sin(i_dec*pi/180.0) ;
     vec_dot = new.L(1)*i_L(1) + new.L(2)*i_L(2) + new.L(3)*i_L(3) ;
     if (vec_dot > 1.0)  
         vec_dot =  1.0 ;
     end
     if (vec_dot < -1.0) 
         vec_dot = -1.0 ;
     end
     arc_dist = acos(vec_dot)*180.0/pi*3600.0 ;
     %/* number = scell2[N_zone-1].star_num[thisnum+i] ; */
     if(abs(arc_dist) < vici && abs(i_mag - new_mag) < mTOL) %/* a possible match */
       if (t > t_start && t < t_end) 
           fprintf(stdout, "Here+-? thisnum=%4d  i=%3d\n", thisnum, i) ;
       end
       thisnum1 = thisnum-i+1 ; 
       for j=0:1:id_count
          if(id_star(j) == stars(scell2(N_zone-1).star_num(thisnum-i)-1).cat_num) 
            astaragain = 1 ;
          end
          if (astaragain ~= 1)  %/* no repeated star */
              intptr  = stars(scell2(N_zone-1).star_num(thisnum-i)-1).cat_num  ;
              fprintf(outdist, "2  %15.6f %5d %15.8f \n", t, intptr, arc_dist) ;
              count = count + 1;
              intptr = intptr + 1;
          end
          i = i + 1;
          astaragain = 0 ;
          if((thisnum - i) < 0)  %/* break ;  */
              thisnum = thisnum + scell2(N_zone-1).num_stars ;  
              if (t > t_start && t < t_end) 
                  fprintf(stdout, "Here3? thisnum=%4d\n", thisnum) ;
              end
          end
          if (t > t_start && t < t_end)
            fprintf(stdout, "+++ ra_diff_BD = %9.2f, ra_diff_new = %9.2f\n",...
                   ra_diff_BD, ra_diff_new) ;
            fprintf(stdout, "  i = %3d  <<<<  i_con = %6.2f\n", ...
               i, scell2(N_zone-1).num_stars/2.0) ;
          end
       end
     end
     kkk=kkk+1;
     if kkk > 500
        i=8;
        ra_diff_BD=BD_limit*3600+1800;
        %keyboard
     end
  end
 thisnum = thisnum0 ;
 if (t > t_start && t < t_end) 
     fprintf(stdout, "Here4? thisnum=%4d i=%3d\n", thisnum, i) ;
 end
 if (thisnum == scell2(N_zone-1).num_stars)
           thisnum = 1 ; %/* was 1 */
 end
 i=0 ;
 while (((abs(ra_diff_BD) < BD_limit*3600+1800.)...
          ||  (abs(ra_diff_new) < FOV_limit*3600+1800.))...  
          &&  i < floor(scell2(N_zone-1).num_stars/2.0))
      i_ra  = stars(scell2(N_zone-1).star_num(thisnum+i)-1).ra  ;
      i_dec = stars(scell2(N_zone-1).star_num(thisnum+i)-1).dec ;
      i_mag = stars(scell2(N_zone-1).star_num(thisnum+i)-1).mag ;

      if (t > t_start && t < t_end)
        fprintf(stdout, "2ND: t=%15.6f N_zone=%3d thisnum=%3d\n", t, N_zone, thisnum);
        fprintf(stdout, "i_ra = %9.2f, i_dec = %9.2f, i_mag = %6.2f\n",...
               i_ra, i_dec, i_mag) ;
        fprintf(stdout, "    star_num = %5d (%2d)\n", ...
           stars(scell2(N_zone-1).star_num(thisnum+i)-1).cat_num, i) ;
      end

      ra_diff_new = abs(new_ra - i_ra) ;
      ra_diff_new = abs(new_ra - i_ra)  ;
      ra_diff_BD  = abs(BD_ra  - i_ra)  ;
      if (ra_diff_new > 180.0)  
          ra_diff_new = 360.0 - ra_diff_new ;
      end
      if (ra_diff_BD  > 180.0)  
          ra_diff_BD  = 360.0 - ra_diff_BD  ;
      end
      ra_diff_new = ra_diff_new * cos(new_dec*pi/180.0) * 3600.0 ;
      ra_diff_BD  = ra_diff_BD * cos(BD_dec*pi/180.0) * 3600.0 ;

      if (abs(new_dec) > 80.0 || abs(BD_dec) > 80.0)   %/* 88 */
         ra_diff_BD  = 100.0 ;
         ra_diff_new = 100.0 ;
      end

      i_L(0) = cos(i_dec*pi/180.0)*cos(i_ra*pi/180.0) ;
      i_L(1) = cos(i_dec*pi/180.0)*sin(i_ra*pi/180.0) ;
      i_L(2) = sin(i_dec*pi/180.0) ;

      vec_dot = L(0)*i_L(0) + L(1)*i_L(1) + L(2)*i_L(2) ;
      if (vec_dot > 1.0)  
          vec_dot =  1.0 ;
      end
      if (vec_dot < -1.0) 
          vec_dot = -1.0 ;
      end
      arc_dist = acos(vec_dot)*180.0/pi*3600.0 ;

      %/* number = scell2[N_zone-1].star_num[thisnum-i] ; */

      if(abs(arc_dist) < vici && abs(i_mag - new_mag) < mTOL) %/* a possible match */
        thisnum1 = thisnum+i+1 ;
        for j=0:1:id_count
           if(id_star(j) == stars(scell2(N_zone-1).star_num(thisnum+i)-1).cat_num) 
             astaragain = 1 ;
           end
        end

        if( astaragain ~= 1)  %/* no repeated star */
          intptr = stars(scell2(N_zone-1).star_num(thisnum+i)-1).cat_num;
          fprintf(outdist, "2  %15.6f %5d %15.8f \n", t, intptr, arc_dist) ;
          count = count+1;
          intptr = intptr + 1;
        end
      end

      i = i + 1;
      astaragain = 0 ;
      if( (thisnum + i) >= scell2(N_zone-1).num_stars) %/* break ; */
          thisnum = thisnum - scell2(N_zone-1).num_stars ;  
          if (t > t_start && t < t_end) 
              fprintf(stdout, "Here5? thisnum=%4d\n", thisnum) ;
          end
      end
 end
 thisnum = thisnum1 ;  

 return;