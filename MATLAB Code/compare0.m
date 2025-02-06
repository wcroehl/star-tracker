function [NEWBLI, ii, jj, esignal] = compare0(n, m, ij, m_star, BLI, NEWBLI)
%Charles

TOL6 =     50./3600.*3.141592653/180.0 ;     % 50 %

tolerance = TOL6 ; % the tol. of difference in interangles between cat. and mesurement

meas_ang =  acos( m_star(n).L(1)*m_star(m).L(1) + ...
                  m_star(n).L(2)*m_star(m).L(2) + ...
                  m_star(n).L(3)*m_star(m).L(3) )*180.0/pi*3600.0; %the cosine of interangle for measurement 
 
for i=1:ij
    for j=(i+1):ij
        cata_ang = acos( BLI(i).L(1)*BLI(j).L(1) + ...
                         BLI(i).L(2)*BLI(j).L(2) + ...
                         BLI(i).L(3)*BLI(j).L(3) )*180.0/pi*3600.0;%the cosine of interangle for catalog 
                                    % of catalog star pairs
        diff = abs(meas_ang - cata_ang) ; % compare the interangles
        
        if(diff < tolerance*180.0/pi*3600.0) %if the angle match 
            [NEWBLI, psignal] = mag_test(n, m, i, j, m_star, BLI, NEWBLI) ;%star mag (britness) match test
            if (psignal <= 2) 
                ii = i;
                jj = j ;
                esignal = psignal; 
                return ;
            end % mag_test was successful(2) or ambiguous(0)
        end
    end
end

ii = i;
jj = i;
%jj = j; %change jj code 11/17
esignal = 3 ;  %/* no match found */
return ;

end
 
