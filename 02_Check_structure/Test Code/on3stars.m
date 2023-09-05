function [NEWBLI, jj, esignal]= on3stars(ij, m_star, BLI, NEWBLI)

TOL3 =     50./3600.*3.141592653/180.0 ;     % 50 %
tolerance = TOL3 ;
dist_tol  = DIST_TOL ;

[jj, sign] = compare2(1, 3, ij, tolerance, m_star, BLI, NEWBLI);% 3rd star search 
if(sign ~= 2)
    [NEWBLI, sign] = compare3(2, 4, jj, tolerance, m_star, BLI, NEWBLI);
end
while(sign == 2  && jj < ij)
    [jj, sign] = compare4(1, 3, jj, ij, m_star, BLI, NEWBLI);
    if(sign ~= 2) 
        [NEWBLI, sign] = compare3(2, 4, jj, tolerance, m_star, BLI, NEWBLI);
    end
end
if(sign ~= 2)
    ambi = sym_test(dist_tol, NEWBLI);
    if (ambi == 0)
        esignal = 7; % 3 stars (all)  
    else
        esignal = 3;
    end
end
if(sign == 2) % 3rd not matched.
    esignal = 3 ;
end
return ;

