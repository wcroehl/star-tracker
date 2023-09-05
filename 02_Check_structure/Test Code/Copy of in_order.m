% int in_order(crf_count, obs_count, Mp, obs_ra, obs_dec, obs_mag, x, y, xmag)
function [c_star, obs_ra, obs_dec, obs_mag] = in_order(crf_count, obs_count, Mp, obs_ra, obs_dec, obs_mag, x, y, xmag, c_star, m_star, crf)
% int    crf_count, obs_count ;
% double Mp(3)(3), *obs_ra, *obs_dec, *obs_mag, x(Ns), y(Ns), xmag(Ns) ;
%  
%  int    i, j, k ;
%  double temporary, *crf_ra, *crf_dec, *crf_mag, tmp_x, tmp_y, tmp_mag ;
%  struct VecStar  *temp  ;
%  struct VecStar2 *temp2 ;

%  temp = (struct VecStar *) malloc(sizeof(struct VecStar)) ; 
%  if ( temp == NULL) 
%    fprintf(stdout, " FAIL : malloc temp (starID_dm:in_order) \n") ; 
%    exit(6) ; 
%  temp2 = (struct VecStar2 *) malloc(sizeof(struct VecStar2)) ; 
%  if ( temp2 == NULL) 
%    fprintf(stdout, " FAIL : malloc temp (starID_dm:in_order) \n") ; 
%    exit(6) ; 
%  crf_ra  = (double *) malloc((crf_count) * sizeof(double) ) ;
%  if ( crf_ra == NULL) 
%    fprintf(stdout, " FAIL : malloc crf_ra (starID_dm:in_order) \n") ; 
%    exit(6) ; 
%  crf_dec = (double *) malloc((crf_count) * sizeof(double) ) ;
%  if ( crf_dec == NULL) 
%    fprintf(stdout, " FAIL : malloc crf_dec (starID_dm:in_order) \n") ; 
%    exit(6) ; 
%  crf_mag = (double *) malloc((crf_count) * sizeof(double) ) ;
%  if ( crf_mag == NULL) 
%    fprintf(stdout, " FAIL : malloc crf_mag (starID_dm:in_order) \n") ;
%    exit(6) ; 

 for(i=1:obs_count)  % compute (RA,DEC) of observed stars %
   for(j=1:3)
      c_star(i).L(j) = 0.0 ;
   end
   for(j=1:3)
       for(k=1:3)
           c_star(i).L(j) = c_star(i).L(j) + Mp(k:j) * m_star(i).L(k) ;
       end
   end
   c_star(i).mag = m_star(i).mag ;

   obs_ra(i) = atan2(c_star(i).L(1),c_star(i).L(0))*180.0/pi ;
                  if ( obs_ra(i) < 0.0 ); obs_ra(i) = obs_ra(i) + 360.0 ; end
   obs_dec(i) = asin(c_star(i).L(2))*180.0/pi ;
   obs_mag(i) = c_star(i).mag ;
 end % the end of for(i<obs_count) %

 for(i=1:obs_count)  
   for(j=i+1:obs_count)  % smaller first %
      if(obs_ra(j) < obs_ra(i))  
        
        temporary = obs_ra(i) ;
        obs_ra(i) = obs_ra(j) ;
        obs_ra(j) = temporary ;

        temporary = obs_dec(i) ;
        obs_dec(i) = obs_dec(j) ;
        obs_dec(j) = temporary ;

        temporary = obs_mag(i) ;    
        obs_mag(i) = obs_mag(j) ;
        obs_mag(j) = temporary ;

        temp2.mag  = c_star(i).mag ; 
        for(k=1:3); temp2.L(k) = c_star(i).L(k) ; end

        c_star(i).mag  = c_star(j).mag ; 
        for(k=1:3); c_star(i).L(k) = c_star(j).L(k) ; end

        c_star(j).mag  = temp2.mag ; 
        for(k=1:3); c_star(j).L(k) = temp2.L(k) ; end

        temp2.mag  = m_star(i).mag ;  
        for(k=1:3); temp2.L(k) = m_star(i).L(k) ; end
	tmp_x = x(i) ;  tmp_y = y(i) ; tmp_mag = xmag(i) ;

        m_star(i).mag  = m_star(j).mag ; 
        for(k=1:3); m_star(i).L(k) = m_star(j).L(k) ; end
	x(i) = x(j) ;   y(i) = y(j) ; xmag(i) = xmag(j) ;

        m_star(j).mag  = temp2.mag ; 
        for(k=1:3); m_star(j).L(k) = temp2.L(k) ; end
	x(j) = tmp_x ;  y(j) = tmp_y ; xmag(j) = tmp_mag ;
      end
   end
 end
    

 for(i=1:crf_count)
   crf_ra(i)  = atan2(crf(i).L(2),crf(i).L(1))*180.0/pi ;
      if ( crf_ra(i) < 0.0 ) crf_ra(i) = crf_ra(i) + 360.0 ; end
   crf_dec(i) = asin(crf(i).L(3))*180.0/pi ;
   crf_mag(i) = crf(i).mag ;
 end
   

 for(i=1:crf_count) 
   for(j=i+1:crf_count)  % smaller IBL first %
      if(crf_ra(j) < crf_ra(i))  
        
        temporary = crf_ra(i) ;
        crf_ra(i) = crf_ra(j) ;
        crf_ra(j) = temporary ;

        temporary = crf_dec(i) ;
        crf_dec(i) = crf_dec(j) ;
        crf_dec(j) = temporary ;

        temporary = crf_mag(i) ; 
        crf_mag(i) = crf_mag(j) ;
        crf_mag(j) = temporary ;

        temp.mag  = crf(i).mag ; 
        for(k=1:3) temp.L(k) = crf(i).L(k) ; end
        temp.num  = crf(i).num ; 

        crf(i).mag  = crf(j).mag ; 
        for(k=1:3) crf(i).L(k) = crf(j).L(k) ; end
        crf(i).num  = crf(j).num ; 

        crf(j).mag  = temp.mag ; 
        for(k=1:3) crf(j).L(k) = temp.L(k) ; end
        crf(j).num  = temp.num ; 
      end
   end
 end
   
 return_value = 1 ;
 