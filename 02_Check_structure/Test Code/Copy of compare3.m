function [NEWBLI, sign] = compare3(n, m, jj, sign, tolerance, m_star, BLI, NEWBLI)
    
meas_ang =  acos( m_star(n).L(1)*m_star(m).L(1) + ...
                  m_star(n).L(2)*m_star(m).L(2) + ...
                  m_star(n).L(3)*m_star(m).L(3) )*180.0/pi*3600.0; %  the cosine of interangle
              
cata_ang = acos( BLI(NEWBLI(2).id).L(1)*BLI(jj).L(1) + ...
                 BLI(NEWBLI(2).id).L(2)*BLI(jj).L(2) + ...
                 BLI(NEWBLI(2).id).L(3)*BLI(jj).L(3) )*180.0/pi*3600.0; 
                 % of catalog star pairs 
                 
diff = abs(meas_ang - cata_ang) ; % compare the interangles 

if(diff < tolerance*180.0/pi*3600.0 )
    NEWBLI(m).L = BLI(jj).L;
    NEWBLI(m).num = BLI(jj).starnum ;
    NEWBLI(m).mag  = BLI(jj).mag ;
    NEWBLI(m).id = jj ;
    return;
end
sign = 2 ;  %/*  not confirmed */
end

%