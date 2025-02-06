function [NEWBLI, psignal] = mag_test(n, m, i, j, m_star, BLI, NEWBLI)
%Charles
%int  n, m, i, j, *psignal ;

Limit =    0.4 ; % 1.0 magnitude limit : starID_pm() %

    meas0 = m_star(n).mag ; %set star magnitude for measurement and catalog
    meas1 = m_star(m).mag ;
    cata0 = BLI(i).mag ;
    cata1 = BLI(j).mag ;
    
    %check whether the two catalog stars i and j must have a magnitude 
    % difference greater than the magnitude error bound ans 
    % if matched star's magnitude is close enough 
    if ( abs(cata0-meas0) < Limit && abs(cata1-meas1) < Limit) 
        if ( abs(cata0-meas1) < Limit && abs(cata1-meas0) < Limit) 
            psignal=0 ; 
        else
            psignal=2 ;
        end
        %save match data to NEWBLI
        NEWBLI(n).L = BLI(i).L;
        NEWBLI(m).L = BLI(j).L;
        NEWBLI(n).num = BLI(i).starnum ;
        NEWBLI(m).num = BLI(j).starnum ;
        NEWBLI(n).mag = BLI(i).mag ;
        NEWBLI(m).mag = BLI(j).mag ;
        NEWBLI(n).id = i ;
        NEWBLI(m).id = j ;
        return ; 
    end

    if ( abs(cata0-meas1) < Limit && abs(cata1-meas0) < Limit)
        if ( abs(cata0-meas0) < Limit && abs(cata1-meas1) < Limit) 
            psignal=0 ;  % practically, this case wouldn't happen %
        else
            psignal=2 ;
        end
        %save match data to NEWBLI
        NEWBLI(n).L = BLI(j).L;
        NEWBLI(m).L = BLI(i).L;
        NEWBLI(n).num = BLI(j).starnum ;
        NEWBLI(m).num = BLI(i).starnum ;
        NEWBLI(n).mag = BLI(j).mag ;
        NEWBLI(m).mag = BLI(i).mag ;
        NEWBLI(n).id = j ;
        NEWBLI(m).id = i ;
        return ; 
    end
    
    psignal = 3 ;   % if magnitude test is failed %
    return ; 
end

