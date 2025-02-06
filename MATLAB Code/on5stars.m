%void on5stars(t, jj, ij, esignal)
function [NEWBLI, jj, esignal] = on5stars(ij, esignal, m_star, BLI, NEWBLI)

TOL6 =     50./3600.*3.141592653/180.0 ;     % 50 %
tolerance = TOL6 ;
% 3rd star searching 
[jj, sign] = compare2(1, 3, ij, tolerance, m_star, BLI, NEWBLI);

if(sign ~= 2)
    [NEWBLI, sign] = compare3(2, 3, jj, sign, tolerance, m_star, BLI, NEWBLI);%if 3rd Match case
end
while(sign == 2  && jj < ij) %if not 3rd match case from compare2 
    [jj, sign] = compare4(1, 3, jj, ij, m_star, BLI, NEWBLI);
    if(sign ~= 2) 
        [NEWBLI, sign] = compare3(2, 3, jj, sign, tolerance, m_star, BLI, NEWBLI);%if 3rd Match case 
    end
end

if (sign ~= 2)   % 3rd star matched, 4th star search 
    [jj, sign] = compare2(1, 4, ij, tolerance, m_star, BLI, NEWBLI);
    if(sign ~= 2) 
        [NEWBLI, sign] = compare3(2, 4, jj, sign, tolerance, m_star, BLI, NEWBLI);%if 4th Match case
    end
    while(sign == 2  && jj < ij) %if not 4th match case from compare2
        [jj, sign] = compare4(1, 4, jj, ij, m_star, BLI, NEWBLI);
        if (sign ~= 2)
            [NEWBLI, sign] = compare3(2, 4, jj, sign, tolerance, m_star, BLI, NEWBLI);%if 4th Match case
        end
    end
    if (sign ~= 2) % 3rd, 4th matched, 5th star search 
        [jj, sign] = compare2(1, 5, ij, tolerance, m_star, BLI, NEWBLI);
        if (sign ~= 2)
            [NEWBLI, sign] = compare3(2, 5, jj, sign, tolerance, m_star, BLI, NEWBLI);            
        end
        while (sign == 2  && jj < ij)
            [jj, sign] = compare4(1, 5, jj, ij, m_star, BLI, NEWBLI);
            if(sign ~= 2)
                [NEWBLI, sign] = compare3(2, 5, jj, sign, tolerance, m_star, BLI, NEWBLI);                
            end
        end
        if (sign ~= 2)
            esignal = 5 ; % 5 stars (all) 
        end
        if (sign == 2)  
            esignal = 4 ;  % 4 stars (0,1,2,3) 
        end
        return ;
    end
    if (sign == 2)   % 3rd matched, 4th not, 5th star search 
        [jj, sign] = compare2(1, 5, ij, tolerance, m_star, BLI, NEWBLI);
        if (sign ~= 2) 
            [NEWBLI, sign] = compare3(2, 5, jj, sign, tolerance, m_star, BLI, NEWBLI);
        end
        while (sign == 2  && jj < ij)
            [jj, sign] = compare4(1, 5, jj, ij, m_star, BLI, NEWBLI);
            if (sign ~= 2) 
                [NEWBLI, sign] = compare3(2, 5, jj, sign, tolerance, m_star, BLI, NEWBLI);
            end
        end
        if (sign ~= 2)  % 3rd OK, 4th not, 5th matched 
            esignal = 4 ;     % 4 stars (0,1,2,4) 
            return ; 
        end
        if (sign == 2)  % 3rd OK, 4th, 5th not matched 
            sign = 3 ;
        end
    end
end

if (sign == 2)   % 3rd star not matched, 4th star search 
    [jj, sign] = compare2(1, 4, ij, tolerance, m_star, BLI, NEWBLI);
    if (sign ~= 2) 
        [NEWBLI, sign] = compare3(2, 4, jj, sign, tolerance, m_star, BLI, NEWBLI);
    end
    while(sign == 2  && jj < ij)
        [jj, sign] = compare4(1, 4, jj, ij, m_star, BLI, NEWBLI);
        if (sign ~= 2)
            [NEWBLI, sign] = compare3(2, 4, jj, sign, tolerance, m_star, BLI, NEWBLI);
        end
    end
    if (sign ~= 2)   % 3rd not matched, 4th matched, 5th search 
        [jj, sign] = compare2(1, 5, ij, tolerance, m_star, BLI, NEWBLI);
        if (sign ~= 2) 
            [NEWBLI, sign] = compare3(2, 5, jj, sign, tolerance, m_star, BLI, NEWBLI);            
        end
        while(sign == 2  && jj < ij)
            [jj, sign] = compare4(1, 5, jj, ij, m_star, BLI, NEWBLI);
            if (sign ~= 2) 
                [NEWBLI, sign] = compare3(2, 5, jj, sign, tolerance, m_star, BLI, NEWBLI);
            end
        end
        if (sign ~= 2)  % 3rd not matched, 4th 5th matched */
            esignal = 4 ;   % 4 stars (0,1,3,4) 
            return ; 
        end
        if (sign == 2)  % 3rd not, 4th OK, 5th not  
            sign = 3 ;
        end
    end
    if (sign == 2)  % 3rd 4th not matched. how about 5th 
        [jj, sign] = compare2(1, 5, ij, tolerance, m_star, BLI, NEWBLI);
        if (sign ~= 2)
            [NEWBLI, sign] = compare3(2, 5, jj, sign, tolerance, m_star, BLI, NEWBLI);            
        end
        while (sign == 2  && jj < ij)
            [jj, sign] = compare4(1, 5, jj, ij, m_star, BLI, NEWBLI);
            if (sign ~= 2) 
                [NEWBLI, sign] = compare3(2, 5, jj, sign, tolerance, m_star, BLI, NEWBLI);
            end
        end
        if (sign ~= 2) % 3rd, 4th not matched, 5th matched 
            sign = 3 ;
        end
    end
end

if (esignal == 0)%If 1st and 2nd star angle match but mag for first and secaond and other stars not match 
     %flip first and secaond star information in NEWBLI
      for k=1:3 
         NBLI.L(k) = NEWBLI(1).L(k) ;
      end
         NBLI.num = NEWBLI(1).num ;
         NBLI.mag   = NEWBLI(1).mag ;
         NBLI.id = NEWBLI(1).id ;

      for k=1:3 
         NEWBLI(1).L(k) = NEWBLI(2).L(k) ;
      end
         NEWBLI(1).num = NEWBLI(2).num ;
         NEWBLI(1).mag   = NEWBLI(2).mag ;
         NEWBLI(1).id = NEWBLI(2).id ;

      for k=1:3 
         NEWBLI(2).L(k) = NBLI.L(k) ;
      end
         NEWBLI(2).num = NBLI.num ;
         NEWBLI(2).mag   = NBLI.mag ;
         NEWBLI(2).id = NBLI.id ;
    
    esignal = esignal + 1;
    return ;
end

%If 1st and 2nd star angle and mag match but other stars not match or other case
if (esignal <= 2)
    esignal = 3 ; 
    return ;
end
