function [jj, sign] = compare2(n, m, ij, tolerance, m_star, BLI, NEWBLI)

Limit =    0.4 ; % 1.0 magnitude limit : starID_pm() %

meas_ang =  acos( m_star(n).L(1)*m_star(m).L(1) + ...
                  m_star(n).L(2)*m_star(m).L(2) + ...
                  m_star(n).L(3)*m_star(m).L(3) )*180.0/pi*3600.0; %  the cosine of interangle

for j=1:ij
    if(    NEWBLI(1).id ~= j && NEWBLI(2).id ~= j ...
        && NEWBLI(3).id ~= j && NEWBLI(4).id ~= j ...
        && NEWBLI(5).id ~= j && NEWBLI(6).id ~= j)
    
        cata_ang = acos( BLI(NEWBLI(1).id).L(1)*BLI(j).L(1) + ...
                         BLI(NEWBLI(1).id).L(2)*BLI(j).L(2) + ...
                         BLI(NEWBLI(1).id).L(3)*BLI(j).L(3) )*180.0/pi*3600.0; 
                         % of catalog star pairs
        
        diff = abs(meas_ang - cata_ang) ; % compare the interangles
        
        if(diff < tolerance*180.0/pi*3600.0)
            mag_diff = abs(BLI(j).mag - m_star(m).mag) ;
            if(mag_diff < Limit)  
                jj = j ;
                sign = 1 ;
                return ;
            end
        end
    end
end
jj = j ;
sign = 2 ;  
end
