function [jj, sign] = compare4(n, m, jj, ij, m_star, BLI, NEWBLI)

Limit =    0.4 ; % 1.0 magnitude limit : starID_pm() %

tolerance = 50./3600.*3.141592653/180.0;
meas_ang =  acos( m_star(n).L(1)*m_star(m).L(1) + ...
                  m_star(n).L(2)*m_star(m).L(2) + ...
                  m_star(n).L(3)*m_star(m).L(3) )*180.0/pi*3600.0; %the cosine of interangle for measurement 
for j=jj+1:ij %loop for diffrence set
    if(    NEWBLI(1).id ~= j && NEWBLI(2).id ~= j ...
        && NEWBLI(3).id ~= j && NEWBLI(4).id ~= j ...
        && NEWBLI(5).id ~= j && NEWBLI(6).id ~= j)  

        cata_ang = acos( BLI(NEWBLI(1).id).L(1)*BLI(j).L(1) + ...
                         BLI(NEWBLI(1).id).L(2)*BLI(j).L(2) + ...
                         BLI(NEWBLI(1).id).L(3)*BLI(j).L(3) )*180.0/pi*3600.0; %catlog star pair angles
        diff = abs(meas_ang - cata_ang) ; % compare the interangles 
        if(diff < tolerance*180.0/pi*3600.0) %If star angles match
            mag_diff = abs(BLI(j).mag - m_star(m).mag) ; %calculate mag diff
            if(mag_diff < Limit) %if mag match
                jj = j ;
                sign = 1 ;
                return ;
            end
        end
    end
end

jj = j ;
sign = 2 ;  % not found 
return ;
