function [t , id_count, m_star, BLI, NEWBLI, return_value, b_star, crf, nstar] = starID_pm(t, q0, nstar, id_count, bore_body, outi, outr, outo, outmag, ixy, x, y, xmag, adjcell, scell, stars, m_star, b_star_prev, crf_prev, rtn_prev,outes1, outes2)
%Charles
% add code



 % Tracker(IST) to Body(IST) here %
  T_B = [ 1.0,    0.0,    0.0
          0.0,    1.0,    0.0
          0.0,    0.0,    1.0 ] ; 

NUMs = 6;
ij = 0 ;
rtn = 0 ;
id_count = 0 ;

Ap = q_to_A(q0); %CHECKED

bore_iner = Ap'*bore_body;

est_del = asin(bore_iner(3))/pi*180.0 ;              
est_alp = atan2(bore_iner(2), bore_iner(1))/pi*180.0;
if(bore_iner(2) < 0.0) est_alp = est_alp + 360.0 ; end

nb_cell = CellIndex(est_alp, est_del) ; % cell num. associated with BD  %CHECKED
                                     % -1 added since array start with zero 
                                     
if (nb_cell == -1)
    error('starID_pm():cell() did not finish normally.') ;
    return_value = 6;
end

BD = cartesian(est_alp*pi/180, est_del*pi/180);  %CHECKED
for i=1:1:adjcell(nb_cell).num_adj
    cellnum = adjcell(nb_cell).cell_num(i)-1; % read adjacent cell # 
    for j=1:1:scell(cellnum).num_stars % read all stars in the cell 
        starnum = scell(cellnum).star_num(j); %catalog 
        al_rad = stars(starnum).ra  / 180 * pi ;  % RA  in radian 
        de_rad = stars(starnum).dec / 180 * pi ;  % Dec in radian 
        L = cartesian(al_rad,de_rad) ;  % get star pos. in cartesian 
        
        BLI(ij+1).IBL = BD'*L;
        
        % dot product of BD and L --> the angle of cosine 
        BLI(ij+1).starnum = starnum+1 ;  % add 1 from array number
        BLI(ij+1).L = L;
        BLI(ij+1).mag  = stars(starnum).mag ; 
        ij = ij+1; % the total # of candidate stars from star catalog */
    end
end

LBLI=length(BLI); %add code from lines 49 to 52, 7/14/23
for i=1:LBLI
    BLI(i).starnum = BLI(i).starnum-1;
end

BLI = sort_BLI(ij, BLI); % structure BLI reorganized with the angular  %CHECKED (starnum_matlab = starnum_c + 1)
                         % distance between BD and star pos. smaller first
nro = nstar - 1;         % the possible # of rotations of measured stars 
k = 0; 
esignal = 3; 
k2 = 0; 
while (esignal == 3 && k < nro)
    for i=1:1:NUMs
        NEWBLI(i).mag = 0.0 ;
        NEWBLI(i).id = -999 ;
        NEWBLI(i).num = 0 ;  
    end
    
    [NEWBLI, ii, jj, esignal] = compare0(1, 2, ij, m_star, BLI, NEWBLI);  %CHECKED
    
    %add fprintf, 07/28/23
    fprintf(outes1, '%15.6f %15.6f %15.6f %15.6f \n',t, ii,jj,esignal);

    i0 = ii ;
    j0 = jj ;
    
    if (esignal == 3)
        [x, y, xmag, m_star] = obs_ang(nstar, x, y, xmag, m_star);
        k = k+1;
    end
end
if (esignal < 3)
    while (esignal < 3) 
        switch(nstar)
            case 6  
                [NEWBLI, jj, esignal] = on6stars(ij, esignal, m_star, BLI, NEWBLI);
                break ;
            case 5  
                [NEWBLI, jj, esignal] = on5stars(ij, esignal, m_star, BLI, NEWBLI);
                break ;
            case 4 
                [NEWBLI, jj, esignal] = on4stars(ij, esignal, m_star, BLI, NEWBLI);  % CHECKED
                  
                %add fprintf, 07/28/23
                fprintf(outes2, '%15.6f %15.6f %15.6f \n',t, jj, esignal);

                break ;
            case 3 
                [NEWBLI, jj, esignal] = on3stars(ij, esignal, m_star, BLI, NEWBLI);
        end
        if (esignal == 1) 
            ii = j0 ;
            jj = i0 ;
        end
        while (esignal == 3 & k < nro && i0 < ij && j0 < ij)
            ii = i0 ;
            jj = j0 ; 
            for i=1:1:NUMs
                NEWBLI(i).mag = 0.0 ;
                NEWBLI(i).id = -999 ;
                NEWBLI(i).num = 0 ;
            end
            [NEWBLI, ii, jj, esignal] = compare1(1, 2, ij, m_star, BLI, NEWBLI);
            
            while (esignal == 3 && k < nro && ii == ij && jj == ij)
                [x, y, xmag, m_star] = obs_ang(nstar, x, y, xmag, m_star);
                k = k+1;
                [NEWBLI, ii, jj, esignal] = compare0(1, 2, ij, m_star, BLI, NEWBLI);
            end
        end
        %  More search for base stars when nstar = 6 
        if (esignal == 3 && k == nro && nstar == 6 && k2 == 0 && ii == ij && jj == ij) 
            [x, y, xmag, m_star] = rearrange(nstar, x, y, xmag, m_star); % Where does it come from?, 7/15/23
            k2=2 ; k=1 ;
        end
        while (esignal == 3 && k < nro && nstar == 6 && k2 > 1 && ii == ij && jj == ij)
            [x, y, xmag, m_star] = obs_ang(nstar, x, y, xmag, m_star);
            k = k+1;
            [NEWBLI, ii, jj, esignal] = compare0(n, m, ij, m_star, BLI, NEWBLI);
            k2 = 5 ;
        end
        if (k2 ~= 0) k2 = -1 ; end
        i0 = ii ;
        j0 = jj ;
    end
end
%esignal
if (esignal == 3) nstar = 0 ; end
if(nstar == 0)
   b_star = b_star_prev;
   crf = crf_prev;
   return_value = rtn_prev
else
    for i=1:nstar
        if (NEWBLI(i).mag ~= 0.0)

            crf(id_count+1).mag   = NEWBLI(i).mag ;
            pos                   = NEWBLI(i).L ;
            crf(id_count+1).num   = NEWBLI(i).num;
            crf(id_count+1).id    = NEWBLI(i).id ;

            norm = NEWBLI(i).L(1)*pos(1)+NEWBLI(i).L(2)*pos(2)+NEWBLI(i).L(3)*pos(3);

            if (norm >  1.0); norm =  1.0 ; end
            if (norm < -1.0); norm = -1.0 ; end
            ang_dist = acos(norm)*180.0/pi*3600.0;
            % angular dist. between before&after aberration correction

            % aberration not simulated 
            crf(id_count+1).L = pos; 

            % for the aberration correction 
            crf(id_count+1).L    = NEWBLI(i).L;

            m_star(id_count+1).L = m_star(i).L;
            x(id_count+1) = x(i);
            y(id_count+1) = y(i);

            % ccd frames to Body-fixed frame 
            b_star(id_count+1).L = zeros(3,1) ;
            for j=1:3
                for k=1:3
                    b_star(id_count+1).L(j) = b_star(id_count+1).L(j) + T_B(j,k) * m_star(id_count+1).L(k);
                end
            end

            norm = sqrt(b_star(id_count+1).L(1)*b_star(id_count+1).L(1)...
                 +      b_star(id_count+1).L(2)*b_star(id_count+1).L(2)...
                 +      b_star(id_count+1).L(3)*b_star(id_count+1).L(3)) ;

            for j=1:3
                b_star(id_count+1).L(j) = b_star(id_count+1).L(j) / norm ; 
            end

            b_star(id_count+1).mag  = m_star(i).mag ; 
            m_star(id_count+1).mag  = m_star(i).mag ; 
            xmag(id_count+1) = xmag(i);

            fprintf(outmag, '%15.3f  1  %6d  %6.2f %6.2f %8.4f\n', ...
                t, crf(id_count+1).num, crf(id_count+1).mag, b_star(id_count+1).mag, ang_dist);

            fprintf(outi, '%15.6f  %20.15f %20.15f %20.15f %10.2f       %10d\n',...
                t, crf(id_count+1).L(1), crf(id_count+1).L(2), crf(id_count+1).L(3), ...
                crf(id_count+1).mag, crf(id_count+1).num); 
           % final way: change to crf(id_count+1).num-1

            fprintf(outo, '%15.6f  %20.15f %20.15f %20.15f %10.2f\n',...
                t, m_star(i).L(1), m_star(i).L(2), m_star(i).L(3), m_star(i).mag);   

            fprintf(ixy,  '%15.6f  %10.2f %10.2f %10.2f %10.2f %10d\n', ...
                t, x(i), y(i), xmag(i), m_star(i).mag, crf(id_count+1).num) ;   

            fprintf(outr, '%15.6f  %20.15f %20.15f %20.15f %10.2f\n', ...
                t, b_star(id_count+1).L(1), b_star(id_count+1).L(2), b_star(id_count+1).L(3),...
                b_star(id_count+1).mag);   
            id_count = id_count + 1;

        end % the end of if (NEWBLI[i].mag != 0.0) 


    end
    return_value = 1;
end 
end

