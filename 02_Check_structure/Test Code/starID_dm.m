function [t, id_count ,ret_fun] = starID_dm(t, q0, nstar, id_count, cctime, w, crf_count, id_star, id_body, outi, outr, outo, outmag, outdist, ixy, x, y, xmag, t_start, t_end)
% FILE   *outi, *outr, *outo,  *outmag, *outdist, *ixy ;
% int    *id_count, nstar,  crf_count, *id_star, *id_body ; 
% double t, q0[4], cctime, w[3], x[Ns], y[Ns], xmag[Ns], t_start, t_end ;
% {  
   %int     i, rtn, in_order(), find(), dmt_write() ; 
   %double  *obs_ra, *obs_dec, *obs_mag, vici, Mp[3][3] ;
   %void    dmt_atti() ;
   if (t > t_start && t < t_end) 
      fprintf("DM: t = %15.6f  q0: %12.8f %12.8f %12.8f %12.8f\n", t, q0(1), q0(2), q0(3), q0(4));
      fprintf("DM: t = %15.6f  w: %12.8f %12.8f %12.8f\n", t, w(1), w(2), w(1));
   end

   Mp = dmt_atti(w, q0);
   
   rtn = in_order(crf_count, nstar, Mp, obs_ra, obs_dec, obs_mag, x, y, xmag) ;
   if (rtn ~= 1)
       fprintf("starID_dm():in_order() didn't finish normally.\n") ;
       ret_fun = 6;
   end
  
   if (t > t_start && t < t_end) 
      fprintf("DM: t = %15.6f  cc_time = %15.6f\n", t, cctime) ;
      fprintf("DM: crf_count = %4d  nstar = %4d  *id_count = %4d\n", ...
		      crf_count, nstar, id_count) ;
   end

   if (t < cctime+5.0)   vici = VICINITY*10.0 ; % /* cctime+1.0,  VI*5 */
   else                  vici = VICINITY+10.0 ; % /* VI+10.0 */ 
   end
   
   if (t > t_start && t < t_end) 
      fprintf("DM: t = %15.6f  cc_time = %15.6f  vici = %8.2f\n",... 
			 t, cctime, vici) ;
   end
   
   rtn = find(t, crf_count, nstar, obs_ra, obs_dec, obs_mag,...
       id_count, id_star, id_body, vici, Mp, outdist, t_start, t_end) ;
   if (rtn ~= 1)
      fprintf("starID_dm():find() didn't finish normally.\n") ;
       ret_fun = 6;
   end

   rtn = dmt_write(t, id_count, id_star, id_body, ...
                outmag, outo, outi, outr, ixy, x, y, xmag) ;
   if (rtn ~= 1)
       fprintf("starID_dm():dmt_write() didn't finish normally.\n") ;
       ret_fun = 6;
   end
   ret_fun = 6;
   return;
% if (t > t_start && t < t_end) 
%       {
%       fprintf(stdout, "DM: t = %15.6f  q0: %12.8f %12.8f %12.8f %12.8f\n", 
% 			 t, q0[0], q0[1], q0[2], q0[3]) ;
%       fprintf(stdout, "DM: t = %15.6f  w: %12.8f %12.8f %12.8f\n", 
% 			 t, w[0], w[1], w[2]) ;
%       }
% 
%    dmt_atti(w, q0, Mp) ;
%    rtn = 0 ;
%    rtn = in_order(crf_count, nstar, Mp, obs_ra, obs_dec, obs_mag, x, y, xmag) ;
%    if (rtn != 1)
%       {
%       fprintf(stdout,"starID_dm():in_order() didn't finish normally.\n") ;
%       exit(6) ; 
%       } 
%   
%    if (t > t_start && t < t_end) 
%       {
%       fprintf(stdout, "DM: t = %15.6f  cc_time = %15.6f\n", 
% 			 t, cctime) ;
%       fprintf(stdout, "DM: crf_count = %4d  nstar = %4d  *id_count = %4d\n", 
% 		      crf_count, nstar, *id_count) ;
%       }
% 
%    if (t < cctime+5.0)   vici = VICINITY*10.0 ;  /* cctime+1.0,  VI*5 */
%    else                  vici = VICINITY+10.0 ;  /* VI+10.0 */ 
%    
%    if (t > t_start && t < t_end) 
%       fprintf(stdout, "DM: t = %15.6f  cc_time = %15.6f  vici = %8.2f\n", 
% 			 t, cctime, vici) ;
%    
%    rtn = find(t, crf_count, nstar, obs_ra, obs_dec, obs_mag, 
%         &(*id_count), id_star, id_body, vici, Mp, outdist, t_start, t_end) ;
%    if (rtn != 1)
%       {
%       fprintf(stdout,"starID_dm():find() didn't finish normally.\n") ;
%       exit(6) ; 
%       } 
% 
%    rtn = dmt_write(t, &(*id_count), id_star, id_body, 
%                 outmag, outo, outi, outr, ixy, x, y, xmag) ;
%    if (rtn != 1)
%       {
%       fprintf(stdout,"starID_dm():dmt_write() didn't finish normally.\n") ;
%       exit(6) ; 
%       } 
% 
%    free ( (void *) obs_mag) ;
%    free ( (void *) obs_ra) ;
%    free ( (void *) obs_dec) ;
% 
%    return(1) ;
