function [NEWBLI, jj, esignal] = on6stars(ij, esignal, m_star, BLI, NEWBLI)
TOL6 =     50./3600.*3.141592653/180.0 ;     % 50 %
tolerance = TOL6 ;
[jj, sign] = compare2(1, 3, ij, tolerance, m_star, BLI, NEWBLI); % 3rd star searching 
if(sign ~= 2); [NEWBLI, sign] = compare3(2, 4, jj, tolerance, m_star, BLI, NEWBLI); end
while(sign == 2  && jj < ij) 
    [jj, sign] = compare4(1, 3, jj, ij, m_star, BLI, NEWBLI);
    if(sign ~= 2); [NEWBLI, sign] = compare3(2, 3, jj, tolerance, m_star, BLI, NEWBLI);   end
end

if (sign ~= 2) % 3rd star matched, 4th star search 
    [jj, sign] = compare2(1, 4, ij, tolerance, m_star, BLI, NEWBLI);
    if(sign ~= 2); [NEWBLI, sign] = compare3(2, 4, jj, tolerance, m_star, BLI, NEWBLI); end
    while(sign == 2  && jj < ij)
        [jj, sign] = compare4(1, 4, jj, ij, m_star, BLI, NEWBLI);
        if (sign ~= 2); [NEWBLI, sign] = compare3(2, 4, jj, tolerance, m_star, BLI, NEWBLI); end
    end
    if (sign ~= 2) % 3rd, 4th matched, 5th star search 
        [jj, sign] = compare2(1, 4, ij, tolerance, m_star, BLI, NEWBLI);
        if (sign ~= 2); [NEWBLI, sign] = compare3(2, 5, jj, tolerance, m_star, BLI, NEWBLI); end
        while (sign == 2  && jj < ij)
            [jj, sign] = compare4(1, 5, jj, ij, m_star, BLI, NEWBLI);
            if(sign ~= 2); [NEWBLI, sign] = compare3(2, 5, jj, tolerance, m_star, BLI, NEWBLI); end
        end
        if (sign ~= 2) % 3rd, 4th, 5th matched, 6th search 
            [jj, sign] = compare2(1, 6, ij, tolerance, m_star, BLI, NEWBLI);
            if (sign ~= 2); [NEWBLI, sign] = compare3(2, 6, jj, tolerance, m_star, BLI, NEWBLI); end
            while (sign == 2  && jj < ij)
                [jj, sign] = compare4(1, 6, jj, ij, m_star, BLI, NEWBLI);
                if(sign ~= 2); [NEWBLI, sign] = compare3(2, 6, jj, tolerance, m_star, BLI, NEWBLI); end
            end
            if (sign ~= 2); esignal = 6 ; end % 6 stars 
            if(sign == 2); esignal = 5 ; end % 5 stars (0,1,2,3,4)
            return; 
        end
        if (sign == 2) % 3rd, 4th matched, 5th not, 6th search */
            [jj, sign] = compare2(1, 6, ij, tolerance, m_star, BLI, NEWBLI);
            if (sign ~= 2); [NEWBLI, sign] = compare3(2, 6, jj, tolerance, m_star, BLI, NEWBLI); end
            while (sign == 2  && jj < ij)
                [jj, sign] = compare4(1, 6, jj, ij, m_star, BLI, NEWBLI);
                if(sign ~= 2); [NEWBLI, sign] = compare3(2, 6, jj, tolerance, m_star, BLI, NEWBLI); end
            end                
            if (sign ~= 2) % 3rd, 4th, 6th matched, 5th not 
                esignal = 5 ;         % 5 stars (0,1,2,3,5) 
                return ; 
            end
            if (sign == 2)  % 3rd, 4th matched, 5th, 6th not */
                esignal = 4 ; % 4 stars (0,1,2,3)   */
                return ; 
            end
        end
    end
    if(sign == 2)  % 3rd matched, 4th not, 5th star search */
        [jj, sign] = compare2(1, 5, ij, tolerance, m_star, BLI, NEWBLI);
        if (sign ~= 2); [NEWBLI, sign] = compare3(2, 5, jj, tolerance, m_star, BLI, NEWBLI); end
        while (sign == 2  && jj < ij)
            [jj, sign] = compare4(1, 5, jj, ij, m_star, BLI, NEWBLI);
            if(sign ~= 2); [NEWBLI, sign] = compare3(2, 5, jj, tolerance, m_star, BLI, NEWBLI); end
        end
        if(sign ~= 2) % 3rd OK, 4th not, 5th matched, 6th search */
            [jj, sign] = compare2(1, 6, ij, tolerance, m_star, BLI, NEWBLI);
            if (sign ~= 2); [NEWBLI, sign] = compare3(2, 6, jj, tolerance, m_star, BLI, NEWBLI); end
            while (sign == 2  && jj < ij)
                [jj, sign] = compare4(1, 6, jj, ij, m_star, BLI, NEWBLI);
                if(sign ~= 2); [NEWBLI, sign] = compare3(2, 6, jj, tolerance, m_star, BLI, NEWBLI); end
            end                
            if (sign ~= 2); esignal = 5; end % 5 stars (0,1,2,3,5) 
            if (sign == 2); esignal = 4; end % 4 stars (0,1,2,3)
            return ; 
        end
        if(sign == 2) % 3rd OK, 4th not, 5th not, 6th search 
            [jj, sign] = compare2(1, 6, ij, tolerance, m_star, BLI, NEWBLI);
            if (sign ~= 2); [NEWBLI, sign] = compare3(2, 6, jj, tolerance, m_star, BLI, NEWBLI); end
            while (sign == 2  && jj < ij)
                [jj, sign] = compare4(1, 6, jj, ij, m_star, BLI, NEWBLI);
                if(sign ~= 2); [NEWBLI, sign] = compare3(2, 6, jj, tolerance, m_star, BLI, NEWBLI); end
            end               
            if (sign ~= 2)  % 3rd OK, 4th not, 5th not, 6th OK
               esignal = 4 ;       % 4 stars (0,1,2,5) 
               return ;
            end
            if (sign == 2)  % 3rd OK, 4th, 5th, 6th not matched 
               % 3 stars among 6 or (undetected) more stars in the FOV
               %   may easily cause ID error */
               sign = 3 ; 
            end
        end
    end
end

if(sign == 2)   % 3rd not matched, 4th star search 
    [jj, sign] = compare2(1, 4, ij, tolerance, m_star, BLI, NEWBLI);
    if(sign ~= 2); [NEWBLI, sign] = compare3(2, 5, jj, tolerance, m_star, BLI, NEWBLI); end
    while(sign == 2  && jj < ij)
        [jj, sign] = compare4(1, 4, jj, ij, m_star, BLI, NEWBLI);
        if (sign ~= 2); [NEWBLI, sign] = compare3(2, 4, jj, tolerance, m_star, BLI, NEWBLI); end
    end
    if (sign ~= 2) % 3rd, 4th matched, 5th star search 
        [jj, sign] = compare2(1, 5, ij, tolerance, m_star, BLI, NEWBLI);
        if (sign ~= 2); [NEWBLI, sign] = compare3(2, 5, jj, tolerance, m_star, BLI, NEWBLI); end
        while (sign == 2  && jj < ij)
            [jj, sign] = compare4(1, 5, jj, ij, m_star, BLI, NEWBLI);
            if (sign ~= 2); [NEWBLI, sign] = compare3(2, 5, jj, tolerance, m_star, BLI, NEWBLI); end
        end
        if (sign ~= 2) % 3rd, 4th, 5th matched, 6th search 
            [jj, sign] = compare2(1, 6, ij, tolerance, m_star, BLI, NEWBLI);
            if (sign ~= 2); [NEWBLI, sign] = compare3(2, 6, jj, tolerance, m_star, BLI, NEWBLI); end
            while (sign == 2  && jj < ij)
                [jj, sign] = compare4(1, 6, jj, ij, m_star, BLI, NEWBLI);
                if(sign ~= 2); [NEWBLI, sign] = compare3(2, 6, jj, tolerance, m_star, BLI, NEWBLI); end
            end
            if (sign ~= 2); esignal = 5 ; end % 6 stars 
            if (sign == 2); esignal = 4 ; end % 5 stars (0,1,2,3,4)
            return; 
        end
        if (sign == 2) % 3rd, 4th matched, 5th not, 6th search */
            [jj, sign] = compare2(1, 6, ij, tolerance, m_star, BLI, NEWBLI);
            if (sign ~= 2); [NEWBLI, sign] = compare3(2, 6, jj, tolerance, m_star, BLI, NEWBLI); end
            while (sign == 2  && jj < ij)
                [jj, sign] = compare4(1, 6, jj, ij, m_star, BLI, NEWBLI);
                if(sign ~= 2); [NEWBLI, sign] = compare3(2, 6, jj, tolerance, m_star, BLI, NEWBLI); end
            end                
            if (sign ~= 2) % 3rd, 4th, 6th matched, 5th not 
                esignal = 4 ;         % 5 stars (0,1,2,3,5) 
                return ; 
            end
            if (sign == 2)  % 3rd, 4th matched, 5th, 6th not */
                esignal = 3 ; % 4 stars (0,1,2,3)   */
                return ; 
            end
        end
    end
    if(sign == 2)  % 3rd matched, 4th not, 5th star search */
        [jj, sign] = compare2(1, 5, ij, tolerance, m_star, BLI, NEWBLI);
        if (sign ~= 2); [NEWBLI, sign] = compare3(2, 5, jj, tolerance, m_star, BLI, NEWBLI); end
        while (sign == 2  && jj < ij)
            [jj, sign] = compare4(1, 5, jj, ij, m_star, BLI, NEWBLI);
            if(sign ~= 2); [NEWBLI, sign] = compare3(2, 5, jj, tolerance, m_star, BLI, NEWBLI); end
        end
        if(sign ~= 2) % 3rd, 4th not, 5th matched, 6th search 
            [jj, sign] = compare2(1, 6, ij, tolerance, m_star, BLI, NEWBLI);
            if (sign ~= 2); [NEWBLI, sign] = compare3(2, 6, jj, tolerance, m_star, BLI, NEWBLI); end
            while (sign == 2  && jj < ij)
                [jj, sign] = compare4(1, 6, jj, ij, m_star, BLI, NEWBLI);
                if(sign ~= 2); [NEWBLI, sign] = compare3(2, 6, jj, tolerance, m_star, BLI, NEWBLI); end
            end                
            if (sign ~= 2); esignal = 4; end % 5 stars (0,1,2,3,5) 
            if (sign == 2); esignal = 3; end % 4 stars (0,1,2,3)
            return ; 
        end
        if(sign == 2) % 3rd, 4th, 5th not matched, how about 6th?
            [jj, sign] = compare2(1, 6, ij, tolerance, m_star, BLI, NEWBLI);
            if (sign ~= 2); [NEWBLI, sign] = compare3(2, 6, jj, tolerance, m_star, BLI, NEWBLI); end
            while (sign == 2  && jj < ij)
                [jj, sign] = compare4(1, 6, jj, ij, m_star, BLI, NEWBLI);
                if(sign ~= 2); [NEWBLI, sign] = compare3(2, 6, jj, tolerance, m_star, BLI, NEWBLI); end
            end               
            if (sign ~= 2); sign = 3; end % 3rd, 4th, 5th not matched, 6th OK
        end
    end
end
      
if (esignal == 0) 
    
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

if (esignal <= 2)
    esignal = 3 ; 
    return ;
end
