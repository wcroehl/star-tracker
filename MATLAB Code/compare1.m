function [NEWBLI, ii, jj, esignal ]= compare1(n, m, ii, jj, ij, m_star, BLI, NEWBLI)

TOL6 =     50./3600.*3.141592653/180.0 ;     % 50 %
tolerance = TOL6 ; % the tol. of difference in interangles between cat. and mesur.

meas_ang =  acos( m_star(n).L(1)*m_star(m).L(1) + ...
                  m_star(n).L(2)*m_star(m).L(2) + ...
                  m_star(n).L(3)*m_star(m).L(3) )*180.0/pi*3600.0; %the cosine of interangle 
              
for i=ii:ij
    if(i > ii) 
        ji = i+1   ;
    else
        ji = jj+1 ;
    end
    if ji > ij
        ji = ij;
        break;
    end
    for j=ji:ij
       cata_ang = acos( BLI(i).L(1)*BLI(j).L(1) + ...
                        BLI(i).L(2)*BLI(j).L(2) + ...
                        BLI(i).L(3)*BLI(j).L(3) )*180.0/pi*3600.0; 
                                    % of catalog star pairs
       diff = abs(meas_ang - cata_ang) ; % compare the interangles */

       if(diff < tolerance*180.0/pi*3600.0 )
           [NEWBLI, psignal] = mag_test(n, m, i, j, m_star, BLI, NEWBLI);      
           if(psignal <= 2)
                ii = i ;
                jj = j ; 
                esignal = psignal ;
                return ; 
           end
       end 
    end
end
ii = i ;   % ij 
jj = j ;   % ij 
esignal = 3 ;
end
