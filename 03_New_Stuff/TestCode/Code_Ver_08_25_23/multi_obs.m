function return_value = multi_obs(t, new, several_num, count, mTOL)
% int  multi_obs(t, new, several_num, count)
% struct VecStar *new ;
% int    *several_num, count ;
% double t ;   % t : not used %
%  {
%  int    i, counted=0, *intptr, iptr, *cptr ; 
%  int    dis_compare() ;
%  double mag_diff ;
counted = 0;

 intptr = 0 ;  % Ns %

 cptr = intptr ;
 for(i=1:count)
   mag_diff=(new.mag-stars((several_num+i)-1).mag) ;
   if (abs(mag_diff) < mTOL)
      counted=counted+1 ;
      cptr = stars((several_num+i)-1).cat_num ;
      cptr=cptr+1 ;
   end
 end

 switch(counted)
   case 0 
       iptr = dis_compare(t, new, several_num, count) ;
       return_value = iptr ;
   case 1 
       iptr = cptr-1 ;
       return_value = iptr ;
   default
   iptr = dis_compare(t, new, intptr, counted) ;
	     return_value = iptr ;
 end

 fprintf(stdout,"starID_dm():find():multi_obs() didn't finish normally.\n") ;
 exit(6) ;
 end