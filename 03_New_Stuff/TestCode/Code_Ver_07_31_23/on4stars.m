function  [NEWBLI, jj, esignal] = on4stars(ij, esignal, m_star, BLI, NEWBLI)

TOL3 =     50./3600.*3.141592653/180.0 ;     % 50 %
tolerance = TOL3 ;

% 3rd star search
[jj, sign] = compare2(1, 3, ij, tolerance, m_star, BLI, NEWBLI);  %CHECKED
if (sign ~= 2) 
    [NEWBLI, sign] = compare3(2, 3, jj, sign, tolerance, m_star, BLI, NEWBLI);    %CHECKED
end
while(sign == 2  && jj < ij) 
    [jj, sign] = compare4(1, 3, jj, ij, m_star, BLI, NEWBLI);  %CHECKING....
    if(sign ~= 2)
        [NEWBLI, sign] = compare3(1, 2, jj, sign, tolerance, m_star, BLI, NEWBLI);
    end
end

if(sign ~= 2)  % 3rd star matched, 4th star search
    [jj, sign] = compare2(1, 4, ij, tolerance, m_star, BLI, NEWBLI);
    if (sign ~= 2)
        [NEWBLI, sign] = compare3(2, 4, jj, sign, tolerance, m_star, BLI, NEWBLI);
    end
    while(sign == 2  && jj < ij)
        [jj, sign] = compare4(1, 4, jj, ij, m_star, BLI, NEWBLI);
        if(sign ~= 2)
            [NEWBLI, sign] = compare3(2, 4, jj, sign, tolerance, m_star, BLI, NEWBLI);
        end
    end
    if(sign ~= 2) % 3rd, 4th star matched 
        esignal = 4 ;   % 4 stars (all) 
        return ; 
    end
    if(sign == 2) % 3rd OK, 4th not 
        sign = 3 ;
    end
end
    
if (sign == 2)  % 3rd star not matched, 4th star search 
    [jj, sign] = compare2(1, 4, ij, tolerance, m_star, BLI, NEWBLI);
    if (sign ~= 2)
        [NEWBLI, sign] = compare3(2, 4, jj, sign, tolerance, m_star, BLI, NEWBLI);
    end
    while(sign == 2  && jj < ij)
        [jj, sign] = compare4(1, 4, jj, ij, m_star, BLI, NEWBLI);
        if(sign ~= 2)
            [NEWBLI, sign] = compare3(2, 4, jj, sign, tolerance, m_star, BLI, NEWBLI);
        end
    end
    if (sign ~= 2) % 3rd not match, 4th star matched */
        sign = 3 ;
    end
end

if (esignal == 0)
    NBLI.L = NEWBLI(1).L;
    NBLI.num = NEWBLI(1).num ;
    NBLI.mag   = NEWBLI(1).mag ;
    NBLI.id = NEWBLI(1).id ;
    
    NEWBLI(1).L = NEWBLI(2).L;
    NEWBLI(1).num = NEWBLI(2).num ;
    NEWBLI(1).mag   = NEWBLI(2).mag ;
    NEWBLI(1).id = NEWBLI(2).id ;
    
    NEWBLI(2).L = NBLI.L;
    NEWBLI(2).num = NBLI.num ;
    NEWBLI(2).mag   = NBLI.mag ;
    NEWBLI(2).id = NBLI.id ;
    
    esignal = esignal + 1;
    return ;
end

if (esignal <= 2)
    esignal = 3 ; 
end
    
    