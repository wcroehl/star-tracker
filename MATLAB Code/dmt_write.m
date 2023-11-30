function [b_star, crf, return_value] = dmt_write(t, id_count, id_star, id_body, outmag, outo, outi, outr, ixy, x, y, xmag, stars, T_B, m_star, b_star, crf)
% int dmt_write(t, id_count, id_star, id_body, outmag, outo, outi, outr, ixy, x, y, xmag)
% FILE   *outmag, *outi, *outr, *outo, *ixy ;
% int    *id_count, *id_star, *id_body  ;
% double t, x(Ns), y(Ns), xmag(Ns) ;
%  {
%  int       i, j, k, rtn=0 ;
%  double    alpha, delta, pos(3), ang_dist ;

for(i=1:id_count)
     
   alpha = stars(id_star(i)).ra  / 180 * pi ;   % RA_in_radian   %
   delta = stars(id_star(i)).dec / 180 * pi ;   % Dec in radian %

   
   pos(1) = cos(delta) * cos(alpha) ; % usually on x axis %
   pos(2) = cos(delta) * sin(alpha) ; % usually on y axis %
   pos(3) = sin(delta) ;  % usually on z axis %
   
   crf(i).L(1) = pos(1) ; % usually on x axis %
   crf(i).L(2) = pos(2) ; % usually on y axis %
   crf(i).L(3) = pos(3) ;  % usually on z axis %
   crf(i).mag  = stars((id_star(i))).mag ;
   crf(i).num  = (id_star(i));

   ang_dist = crf(i).L(1)*pos(1)+crf(i).L(2)*pos(2)+crf(i).L(3)*pos(3) ;
   if (ang_dist >  1.0);  ang_dist =  1.0 ; end
   if (ang_dist < -1.0);  ang_dist = -1.0 ; end
   ang_dist = acos(ang_dist)*180.0/pi*3600.0; 

   for(j=1:3); crf(i).L(j) = pos(j) ; end

   for(j=1:3); b_star(i).L(j) = 0.0 ; end
   for(j=1:3)
       for(k=1:3)
         b_star(i).L(j) = b_star(i).L(j) + T_B(j,k) * m_star(id_body(i)).L(k) ;
       end
   end
   b_star(i).mag = m_star(id_body(i)).mag;

   fprintf(outmag," %15.3f  2  %6d  %6.2f %6.2f %8.4f\n", t, crf(i).num, crf(i).mag, b_star(i).mag, ang_dist) ;

   fprintf(outi,"%15.6f  %20.15f %20.15f %20.15f %10.2f       %10d\n", t, crf(i).L(1), crf(i).L(2), crf(i).L(3), crf(i).mag, crf(i).num) ;
   
   fprintf(outo,"%15.6f  %20.15f %20.15f %20.15f %10.2f\n", t, m_star((id_body(i)+1)).L(1), m_star((id_body(i)+1)).L(2),m_star((id_body(i)+1)).L(3), m_star((id_body(i)+1)).mag) ;
   
   fprintf(ixy, "%15.6f  %10.2f %10.2f %10.2f %10.2f  %10d\n", t, x((id_body(i)+1)), y((id_body(i)+1)), xmag((id_body(i)+1)), m_star((id_body(i)+1)).mag, crf(i).num) ;

   fprintf(outr,"%15.6f  %20.15f %20.15f %20.15f %10.2f\n", t, b_star(i).L(1), b_star(i).L(2), b_star(i).L(3), b_star(i).mag) ;
 end
       
 return_value = 1 ;
 
end