function return_value = dis_compare(t, new, intptr, counted, stars)
% int dis_compare(t, new, intptr, counted)
% struct VecStar *new ;
% int    counted, *intptr ;
% double t ;  % t : not used %
%  {
%  int     i, smallest ;
%  double  i_L(3), vec_dot, arc_dist, arc_dist_s ;

 i_L(0) = cos(stars(intptr-1).dec*pi/180.0)*cos(stars(intptr-1).ra*pi/180.0) ;
 i_L(1) = cos(stars(intptr-1).dec*pi/180.0)*sin(stars(intptr-1).ra*pi/180.0) ;
 i_L(2) = sin(stars(intptr-1).dec*pi/180.0) ;
 vec_dot = new.L(0)*i_L(0) + new.L(1)*i_L(1) + new.L(2)*i_L(2) ;
 if (vec_dot > 1.0);  vec_dot =  1.0 ; end
 if (vec_dot < -1.0); vec_dot = -1.0 ; end
 arc_dist_s = acos(vec_dot) ;
 smallest = intptr ;
     
 for(i=1:counted-1)
     
   i_L(0) = cos(stars((intptr+i)-1).dec*pi/180.0)*cos(stars((intptr+i)-1).ra*pi/180.0) ;
   i_L(1) = cos(stars((intptr+i)-1).dec*pi/180.0)*sin(stars((intptr+i)-1).ra*pi/180.0) ;
   i_L(2) = sin(stars((intptr+i)-1).dec*pi/180.0) ;

   vec_dot = new.L(0)*i_L(0) + new.L(1)*i_L(1) + new.L(2)*i_L(2) ;
   if (vec_dot > 1.0);  vec_dot =  1.0 ; end
   if (vec_dot < -1.0); vec_dot = -1.0 ; end
   arc_dist = acos(vec_dot) ;

   if( abs(arc_dist) < abs(arc_dist_s))
      smallest = intptr+i ;
      arc_dist_s = arc_dist ;
   end
 end
 return_value = smallest ; 
end