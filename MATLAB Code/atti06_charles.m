clear, clc

load('atti06_h.mat');

format long
DEL_TIME = 1000.0;
FNLEN =    256;

now = clock;
tic

%%% Lee
znum = 0;
zsum = [0 0 0];

b_est_sum = [0,0,0];
b_est_cnt = 0;

% INPUT : Gyroscope 
gyro = fopen("D:\StarTracker-20230706T212348Z-001\StarTracker\02_Check_structure\data\g_rate_lrs.dat", 'r'); 

% INPUT : LRS STAR 
sttr = fopen("D:\StarTracker-20230706T212348Z-001\StarTracker\02_Check_structure\data\FOVs.dat", 'r'); 

% INPUT : INITIAL ATTITUDE 

msg = "this is dumb";
% OUTPUTS 
outi = fopen('inertials.dat', 'w'); 
outr = fopen('rotatings.dat', 'w'); 
rate = fopen('gyrorate.dat', 'w');
outM = fopen('ckg_M.dat', 'w');
out1a = fopen('ckg_1a.dat', 'w');
outpos = fopen('glas_pos.dat', 'w');
out1b = fopen('ckg_1b.dat', 'w');
outo = fopen('observed.dat', 'w');
out2 = fopen('ckg_2.dat', 'w');
outds = fopen('rmsstar.dat', 'w');
outdq = fopen('rmsckg.dat', 'w');
outn = fopen('IDed_num.dat', 'w');
out3 = fopen('ckg_3.dat', 'w');
gbias = fopen('ckg_gbias.dat', 'w');
outgbias = fopen('gbias.average', 'w');
angvel = fopen('ckg_angvel.dat', 'w');
outmag = fopen('star_mag.dat', 'w');
outdist = fopen('dm_dist.dat', 'w');
qc = fopen('ist_qc.dat', 'w');
merge = fopen('merge.dat', 'w');
ixy = fopen('ided_xy.dat', 'w');
t_diff = fopen('t_diff.dat', 'w');

% Test Output, added, 07/28/23

% P_theta Test Output
outPT1 = fopen('outpt1.dat', 'w');
outPT2 = fopen('outpt2.dat', 'w');
outPT3 = fopen('outpt3.dat', 'w');

% esignal Test Output
outes1 = fopen('outes1.dat', 'w');
outes2 = fopen('outes2.dat', 'w');

atti_rms=0; %added, 8/24/23
%*****************************************************************/

% load all data in datafile in memory area */
load('star_data.mat')

id_star = zeros(Ns,1);
id_body = zeros(Ns,1);




bore_body = T_B(:,3);  %(0,0,1) in the ST frame 

ccdstep = T_STEP_CCD;
rmssum  = zeros(3,1);
rmslsum = zeros(3,1);
b_add   = zeros(3,1);
n_add = 0 ;
n_g = 0 ;
nsig = 0 ;
rg = 0 ;
consec = 0;
crf_count = 0 ;
cnt = 0 ;
cnt0 = 0 ;
ID_loca = 'F' ;
id_ = 99 ; % Failed as initial: Used actually? 
total_mstar = 0.0 ;
total_ided = 0.0 ;
sign1 = 0;
QC1 = 0;
M = 0;

g_angle = zeros(4,1);
gt_old = 0.0 ;

% KF (t < KF) : in opt. bench coordinate frame 

w0 = [2.0*pi/T/10.0 + 1e-4, 1.e-4, 1.e-4];

inival = dlmread('q_initial.txt')
cit = inival(1,1);
q = inival(2,:);
b_est_average = (0.01/3600*pi/180)*ones(1,3); 
q = q/norm(q);
q0 = q;
qq_old = q;
qq = zeros(4,1);
P_theta = zeros(3);

P = diag([ 1e-6 1e-6 1e-6 1e-8 1e-8 1e-8 ]);
b_star_prev = 0;
crf_prev = 0;
rtn_prev = 0;

% READ STAR DATA 
ccd_time = -2000.0; 
load('measurement_data.mat');  % Data file containing FOVs_meas, gyro_meas, q0, t0
nstar = 0;
FOVs_count = 1;
while (nstar < 3)
    [ccd_time, x, y, xmag, m_star, nstar] = read_CCD0(FOVs_meas, Ns, FOVs_count);  
    FOVs_count = FOVs_count + 1;
end

% READ GYRO DATA %
gyro_count = 1;
[t_gyro, w, u, gread] = read_gyro(b_est_average, gyro_meas, gyro_count);
gyro_count = gyro_count +1;
if gread == 1
    fprintf(' No more gyro data\n');
end

while ( (t_gyro - ccd_time) > 0.05)
    [ccd_time, x, y, xmag, m_star, nstar] = read_CCD(ccd_time,FOVs_meas, Ns);
    if nstar > Ns
        fprintf('Too many stars at t=%8.2f\n', ccd_time);
    elseif nstar == -999
        fprintf('Empty IST data file before "while" loop\n') ;
    else
        total_mstar = total_mstar + nstar ;
    end
    FOVs_count = FOVs_count + 1;
end

while ( (ccd_time - t_gyro) > 0.05) 
    [t_gyro, w, u, gread] = read_gyro(b_est_average, gyro_meas, gyro_count);
    if gread == 1
        fprintf(' Empty gyro data file at ccd_time = %8.2f\n', ccd_time) ;
    end
    gyro_count = gyro_count +1;
end

ccd_time0 = ccd_time ;
cctime = ccd_time + Ctime ; % set up the pm range 
t_rms = t_gyro ;
if ccd_time < KF
    w = w0;
end

c_dm = 0;
c_dm_ided = 0;
c_pm = 0;
c_pm_ided = 0;

while( ccd_time <= TIMELIMIT ) 
    d_i = 0.0 ;    
    fprintf(t_diff, '%15.6f  %15.6f  %15.6f\n', ccd_time, t_gyro, ccd_time-t_gyro) ;
    fprintf('%15.6f  %15.6f  %15.6f\n', ccd_time, t_gyro, ccd_time-t_gyro) ;
    q0 = q0/norm(q0);
    cnt = 0;
    if ( ccd_time < cctime && nstar >= 3)
        c_pm = c_pm + 1 ;
        rtn = 0 ;
        if (nstar > 6)
            nstar = 6 ;
        end % up to 6 stars for PM */
        % rtn = starID_pm(ccd_time, q0, nstar, &cnt, bore_body, outi, outr, outo, outmag, ixy, x, y, xmag) ;
        [ccd_time, cnt ,m_star, BLI, NEWBLI, rtn, b_star, crf, nstar] = starID_pm(ccd_time, q0, nstar, cnt, bore_body, outi, outr, outo, outmag, ixy, x, y, xmag, adjcell, scell, stars, m_star, b_star_prev, crf_prev,rtn_prev, outes1, outes2);   
        %add outes's, 07/28/23
        b_star_prev = b_star;
        crf_prev = crf;
        rtn_prev = rtn;
        fprintf("NEWBLI\n") ;        %  !!PRINT STATEMENT!!<-----------------------------------------------------------------*/
        fprintf("mag        id      num\n") ;
        NEWBLI_size = size(struct2table(NEWBLI),1);
        for i=1:NEWBLI_size
            fprintf("%f  %d   %4d\n",NEWBLI(i).mag,NEWBLI(i).id,NEWBLI(i).num) ;
        end
        fprintf("\n");
        
        fprintf("BLI\n");        %  !!PRINT STATEMENT!!<-----------------------------------------------------------------*/
        fprintf("IBL     starnum    L                             mag\n") ;
        BLI_size = size(struct2table(BLI),1);
        for i=1:BLI_size
            fprintf("%0.4f  %4d      %0.4f   %0.4f   %0.4f   %4f\n",BLI(i).IBL,BLI(i).starnum,BLI(i).L(1),BLI(i).L(2),BLI(i).L(3), BLI(i).mag) ;
        end
        fprintf("\n");
        
        fprintf("m_star\n") ;        %  !!PRINT STATEMENT!!<-----------------------------------------------------------------*/
        fprintf("L                             mag\n") ;
        m_star_size = size(struct2table(m_star),1);
        for i=1:m_star_size
            fprintf("%0.4f   %0.4f   %0.4f   %4f\n",m_star(i).L(1),m_star(i).L(2),m_star(i).L(3),m_star(i).mag) ;
        end
        fprintf("\n");
        
        for i=1:1:Ns 
            x(i) = 0.0;
            y(i) = 0.0;
            xmag(i) = 0.0 ; 
        end
        if rtn ~= 1
            fprintf('starID_pm() finished abnormally at %10.2f\n',ccd_time);
        end
        if (cnt > 0)
            c_pm_ided = c_pm_ided + 1;
        end
        crf_count = cnt ;
        ID_loca = 'P' ;
        id_ = 12 ;
        if (cnt >= 3 && rg > 0 )  % UPDATE from > 0 
            cnt0 = 0;
            n_add = 0;
            n_g = 0;
            rg = 0 ;
            b_add = zeros(3,1);
        end
    elseif (ccd_time >= cctime && nstar > 0 && consec ~= 0)
        c_dm = c_dm + 1;
        [ccd_time, cnt ,rtn] = starID_dm(ccd_time, q0, nstar, cnt, cctime, w, crf_count, id_star, id_body, outi, outr, outo, outmag, outdist, ixy, x, y, xmag, ccd_time0, ccd_time, m_star, crf, stars, T_B, b_star, scell2);
        for i=1:BLI_size
            fprintf("%0.4f  %4d      %0.4f   %0.4f   %0.4f   %4f\n",BLI(i).IBL,BLI(i).starnum,BLI(i).L(1),BLI(i).L(2),BLI(i).L(3), BLI(i).mag) ;
        end
        fprintf("\n");
        
        fprintf("m_star\n") ;        %  !!PRINT STATEMENT!!<-----------------------------------------------------------------*/
        fprintf("L                             mag\n") ;
        m_star_size = size(struct2table(m_star),1);
        for i=1:m_star_size
            fprintf("%0.4f   %0.4f   %0.4f   %4f\n",m_star(i).L(1),m_star(i).L(2),m_star(i).L(3),m_star(i).mag);
        end
        fprintf("\n");
        for i=1:1:Ns 
            x(i) = 0.0;
            y(i) = 0.0;
            xmag(i) = 0.0 ; 
        end
        if rtn ~= 1
            fprintf('starID_dm() finished abnormally at %10.2f\n', ccd_time);
        end
        if (cnt > 0)
            c_dm_ided = c_dm_ided + 1;
        end
        crf_count = cnt ;
        ID_loca = 'D' ;
        id_ = 13 ;
        if (cnt > 1) % UPDATE from 0 
            cnt0 = 0 ;
            if (rg > 0)
                n_add = 0;
                n_g = 0;
                rg = 0 ;
                b_add = zeros(3,1);
            end
        end
    else
        cnt = 0;
    end
    fprintf('%g %g\n', ccd_time, TIMELIMIT); %FOR DEBUGGING
    total_ided = total_ided + crf_count ;  % just for statistics %
 
    
    if (cnt >= 3)
        for i = 1:numel(crf)
            for j = 1:numel(crf(i).L)
                V(i,j) = crf(i).L(j);
                V(i,j+1) = crf(i).mag;
                
                W(i,j) = b_star(i).L(j);
                W(i,j+1) = b_star(i).mag;
            end
        end
        
        %add fprintf, 07/28/23
        fprintf(outPT1, '%15.6f \n%.15f %.15f %.15f \n%.15f %.15f %.15f \n%.15f %.15f %.15f \n',ccd_time, P_theta);

        [qq, P_theta, d_i, rtn] = q_method(W, V, nstar, d_i, outPT2, outPT3, ccd_time);% added outPT's and ccd_time, 07/28/23

        %qq = quest(W(1:3,:),V(1:3,:), ones(nstar,1)/nstar)
        
        rtn = 0;
        if (rtn ~= 0)
			  fprintf(stdout,"q_method() finished abnormally at %10.2f\n", ccd_time);
        end
        % Discontinuity Check
        for i = 1:1:3
            if(qq(i) > .3 && qq(i)*qq_old(i)<0) 
                opp_sign=opp_sign+1;
                if(opp_sign>=2) 
                    sign1=1;
                else
                    sign1=0;
                end
            end
        end  
        opp_sign=0;
        if (sign1 == 1)
            if nsig == 0
                nsig = 1;
            else
                nsig = 0;
            end
            sign1 = 0;
        end
        qq_old = qq;
        if nsig == 1
            qq = -qq;
        end
        
        fprintf(out1a,"%15.6f", ccd_time);
        fprintf(out1a,"  %15.10f  %15.10f  %15.10f  %15.10f",qq);
        fprintf(out1a,"\n");
        fprintf(out2,"%15.6f", ccd_time);
        fprintf(out2,"\t %12.8f \t %12.8f \t %12.8f", sqrt(4*P_theta(1,1))*180./pi*3600,...
              sqrt(4*P_theta(2,2))*180./pi*3600, sqrt(4*P_theta(3,3))*180./pi*3600); %edit code from lines 300 to 301, 7/27/23
       % original code: fprintf(out2,"\t12.8f\t%12.8\t%12.8 ", sqrt(4*P)); 
       % completely diffrent from C, 7/10/23
        fprintf(out2,"\n");
       
       % there are errors? 7/28/23 
        if(consec < 9)
            consec = consec+1;
            if(consec == 9)
               %10th consecutive data: (re)start filter

               %changed, 8/24/23
                P = zeros(6,6);
                for i=1:1:3
                    for j=1:1:3
                        P(i,j) = P_theta(i,j);
                    end
                end

%                for(i=0;i<3;
%                   for(j=0;j<3;j++) P[i][j] = P_theta[i][j] ;

               % P = P_theta; original code
                t_old = ccd_time;
            end
            A = q_to_A(qq);
            qp = q;
            q = q0;
            q0 = qq;    
% 			  for(i=0;i<3;i++) for(j=0;j<3;j++) M[i][j] = 0.0 ;    
            M = zeros(3,3);
% 			  for(i=0;i<3;i++) for(j=0;j<3;j++) for(k=0;k<3;k++)
% 				  M[i][j] += T_B[k][i] * A[k][j]   ;
            M = M + T_B*A;
            fprintf(out1b,"%15.6f 11 ", ccd_time);
            fprintf(out1b,"\n");
            fprintf(out1b,"\t%15.12f",qq);
            fprintf(out1b,"\n");
            fprintf(qc,"%15.6f\t%3d\t%6.2f\n",ccd_time, 1, -0.5);
            fprintf(outn,"%15.6f  11  %2d  %2d  %2d  %8.3f\n", ccd_time, nstar, cnt, consec, 0.0);
            fprintf(outM,"%12.6f %3d\n%20.15f %20.15f %20.15f \n%20.15f %20.15f %20.15f \n% 20.15f %20.15f %20.15f \n",...
                  ccd_time, QC1, inv(M)) ;%updated 7/27 from M
        end
    end
    if (consec == 10 && cnt ~= 0)
        t_elapse = ccd_time - ccd_time0 ;    % Tuning adjustment
		b_est = b_est_average;
% 		  cuvfout = cuvf(ccd_time, t_old, t_elapse, qq, q, P, w, b_est, cnt, zsum, znum, ccdstep); I dont like this
            % errors related to P happen below, 7/28/23
        [zsum, znum, b_est, q, t_old, P] = cuvf(ccd_time, t_old, t_elapse, qq, q, P, w, b_est, cnt, zsum, znum, ccdstep, crf, b_star); %original code
            %[zsum, znum, b_est, q, t_old, P] = cuvf(ccd_time, t_old, t_elapse, qq, q, P, w, b_est, cnt, zsum, znum, ccdstep, crf); %add crf, 08/25/23
        cuvfout = 0;
		b_est_average = b_est;
        if ccd_time >= 0.0                   %  not for me
		    b_est_sum = b_est_sum + b_est;
			b_est_cnt = b_est_cnt+1;
        end
        if cuvfout ~= 0
		    fprintf(stdout," cuvf() finished abnormally at %10.2f\n", ccd_time) ;
			error(msg) ;
        end
        ready_merge = 0;
        if ready_merge == 1
	        fprintf(merge," %15.6f  %15.6f\n", ccd_time-dt_merge/10.0, ccd_time-1.0) ;
	        ready_merge = 0 ;
        end
        if ccd_time < KF  
            w0 = w;
        end
      
		qp = q0; 
        q0 = q;
		A = q_to_A(qp) ;

		M = zeros(3,3);
        for i=1:1:3
            for j=1:1:3 
                for k=1:1:3
                    M(i,j) = T_B(k,i) * A(k,j) ;
                end
            end
        end
        if (ccd_time > cctime)  % for QC          %dont know 5/4/22
            rmssum = rmssum + zsum;
            rmsnum  = rmsnum + znum ;
            rmslsum = rmslsum + zsum;
            rmslnum = rmslnum + znum;
        end
        if (cnt >= 3)
            qp_in = -qp;
            qp_in(3) = qp(3);
            z_bar = qcomp(qq,qp_in,z_bar) ;
            temp = sqrt(z_bar(1)*z_bar(1) + z_bar(2)*z_bar(2)+z_bar(3)*z_bar(3) + z_bar(4)*z_bar(4)) ;
            z_bar = z_bar/temp ;
            z = 2*z_bar; % 3 x 1 vector
            if (t_gyro > cctime)   % for QC
                rmszsum = rmszsum + sqrt(z*z) ;
                rmszscalar = rmszscalar + sqrt(z_bar(4)*z_bar(4));
                rmsznum = rmsznum + 1;
                rmszlsum = rmszlsum + sqrt(z*z) ;
                rmszlscalar = rmszlscalar + sqrt(z_bar(4)*z_bar(4)) ;
                rmszlnum = rmszlnum + 1;
            end
        end
        if ( cnt == 2 ) % || (ID_loca == 'D' && cnt >= 3) )
		    d_i = b_star(1).L(1)*b_star(2).L(1) + b_star(1).L(2)*b_star(2).L(2) + b_star(1).L(3)*b_star(2).L(3);
        end
        if ( cnt  < 2 )
            d_i = 0.0 ;
        end

		fprintf(gbias, "%15.6f  %15.8e  %15.8e  %15.8e\n", ccd_time, b_est*3600/3.1415*180);
		fprintf(angvel, "%15.6f  %15.5f  %15.5f  %15.5f\n", ccd_time, w*180/3.1415*36000);
		OneSig1 = sqrt(P(1,1))*180./pi()*3600;
		OneSig2 = sqrt(P(2,2))*180./pi()*3600;
		OneSig3 = sqrt(P(3,3))*180./pi()*3600;
		fprintf(out3,"%15.6f",ccd_time) ;
		fprintf(out3,"   %12.8f   %12.8f   %12.8f    %12.8f  %12.8f  %12.8f\n", OneSig1, OneSig2, OneSig3, sqrt(P(4,4))*180/pi*3600,sqrt(P(5,5))*180/pi*3600,sqrt(P(6,6))*180/pi*3600);
          %nope
		atti_rms = sqrt(OneSig1*OneSig1+OneSig2*OneSig2) ;

        if atti_rms <= 2.0
            QC1 = 2;
        elseif atti_rms > 2.0 && atti_rms < 4.0
            QC1 = 1;                
        else
            QC1 = 0 ;
        end

        if (ID_loca == 'P')  
            fprintf(out1b,"%15.6f  12 ",ccd_time) ;
        elseif (ID_loca == 'D')  
            fprintf(out1b,"%15.6f  13 ",ccd_time) ;
        end
		fprintf(out1b,"  %15.12f",qp);
		fprintf(out1b,"\n") ;
		fprintf(qc,"%15.6f  %3d %6.2f\n",ccd_time, QC1, atti_rms) ;
		fprintf(outM,"%12.6f %3d\n%20.15f %20.15f %20.15f \n%20.15f %20.15f %20.15f \n% 20.15f %20.15f %20.15f \n", ccd_time, QC1, inv(M));
           %updated 7/27 from M
		
        if (cnt > 2)
		    fprintf(outn,"%15.6f  %2d  %2d  %2d  %2d  %8.3f\n", ccd_time, id_, nstar, crf_count, consec, acos(d_i/cnt)*180.0/pi());
        elseif (cnt == 2)
		    fprintf(outn,"%15.6f  %2d  %2d  %2d  %2d  %8.3f\n", ccd_time, id_, nstar, crf_count, consec, acos(d_i)*180.0/pi());
		else % if (cnt == 1)  */
            fprintf(outn,"%15.6f  %2d  %2d  %2d  %2d  %8.3f\n",ccd_time, id_, nstar, crf_count, consec, 0.0);  
        end  
    end %end of if ( consec == 10 && cnt != 0) */

    if (cnt == 0)
	    cnt0 = cnt0 + 1;
        if (cnt0 > 100.0 || (cnt0 > 10.0 && nstar >= 3) )    % Changed
		    cctime = ccd_time + Ctime ;
			consec = 0 ;
			QC1 = 0 ;
            atti_rms = -1.0 ;
        end
        if (cnt0 > 600.0)
		    % fprintf(stdout, " 1. Warning: No obs for %6.2f min at %15.6f\n",cnt0/10.0/60.0, t_gyro) ; */
			dt_merge = cnt0 ;
			ready_merge = 1 ;
			QC1 = 0 ; atti_rms = -2.0 ;
        end
          
		[rg, q0] = no_stars(q0, w) ; % no idea % make no_stars work, 7/14/23 
		qp = q0;
		A = q_to_A(qp) ;
		M = zeros(3,3);
        for i=1:1:3 
            for j=1:1:3 
                for k=1:1:3
                    M(i,j) = M(i,j) + T_B(k,i) * A(k,j);
                end
            end
        end
		fprintf(out1b,"%15.6f  21 ", ccd_time) ; % no star IDed */
		fprintf(out1b,"  %15.12f",qp) ;
		fprintf(out1b,"\n") ;
		fprintf(qc,"%15.6f  %3d %6.2f\n",ccd_time, QC1, atti_rms) ; % change the code, 8/24/23 %dont have 
		ID_loca = 'G' ; id_ = 21 ; consec = 0 ;
		fprintf(outn,"%15.6f  %2d  %2d  %2d  %2d  %8.3f\n", ccd_time, id_, nstar, cnt, consec, 0.0) ; 
		fprintf(outM,"%12.6f %3d\n%20.15f %20.15f %20.15f \n%20.15f %20.15f %20.15f \n% 20.15f %20.15f %20.15f \n",ccd_time, QC1, inv(M)) ;
           %updated 7/27 from M
    end
	if (consec == 9 && t_old == ccd_time) 
        consec = 10;
    end

    [ccd_time, x, y, xmag, m_star, nstar] = read_CCD(ccd_time,FOVs_meas, Ns);
    while nstar == 0
        if ( nstar == -999 )
		    fprintf(stdout,"No more IST data at or after %10.2f\n",ccd_time);
            if ( t_elapse < DEL_TIME )    
                error(msg) ;
            else
                break;
            end
        end
		[t_gyro, w, u, gread] = read_gyro(t_gyro, u, b_est_average, w, gyro);

        if (gread == 1)
		    fprintf(stdout,"No more gyro data. IST data may exist at %10.2f\n", t_gyro);
            if (t_gyro < 86400.0)  
                error(msg) ;
            else
			    fprintf(outgbias, "%15ld  %15.5f  %15.5f  %15.5f\n", b_est_cnt, b_est_sum/b_est_cnt*3600/PI*180.0) ;
				error(msg) ;
            end
        end
        if (ccd_time < KF)  
            w = w0;
        end
		[rg,q0] = no_stars(q0, w) ; % change 7/14/23 
        %rg = no_stars(ccd_time, t_gyro, q0, w) ;% original code    %dont have
		qp = q0;

		q_to_A(qp, A) ;

		M = zeros(3,3);
        for i=1:1:3 
            for j=1:1:3 
                for k=1:1:3
                    M(i,j) = M(i,j) + T_B(k,i) * A(k,j);
                end
            end
        end
		cnt0 = cnt0 + 1;
        if (cnt0 > 100.0 || (cnt0 > 10.0 && nstar >= 3) ) % Changed */
		    cctime = ccd_time + Ctime ;
			consec = 0 ;
			QC1 = 0 ; atti_rms = -1.0 ;
        end
        if (cnt0 > 600.0)
		    dt_merge = cnt0 ;
			ready_merge = 1 ;
			QC1 = 0 ; atti_rms = -2.0 ;
        end
		fprintf(out1b,"%15.6f  22 ",ccd_time) ; % No star observed */
		fprintf(out1b,"  %15.12f",qp) ;
		fprintf(out1b,"\n") ;

		ID_loca = '2' ;  
        id_ = 22 ; 
        consec = 0 ;
		fprintf(qc,"%15.6f  %3d %6.2f\n",ccd_time, QC1, atti_rms) ;
		fprintf(outn,"%15.6f  %2d  %2d  %2d  %2d  %8.3f\n",ccd_time, id_, 0, 0, consec, 0.0) ;
		fprintf(outM,"%12.6f %3d\n%20.15f %20.15f %20.15f \n%20.15f %20.15f %20.15f \n% 20.15f %20.15f %20.15f \n",ccd_time, QC1, inv(M)) ;
           %updated 7/27 from M
    end % end of while (nstar = read_CCD(  ) ) */

    if (nstar > 0) 
        total_mstar = total_mstar + nstar ;
    end
    if (nstar > Ns )
	    fprintf(stdout,"Too many stars at t=%8.2f\n", ccd_time) ;
		error(msg) ;
    end
    if (nstar == -999) 
        break; % 2nd break : 1st break occurred with no CCD data */
    end
 	[t_gyro, w, u, gread] = read_gyro(b_est_average, gyro_meas, gyro_count);
    if (gread == 1)
 	    fprintf(stdout, "No more gyro data at %10.2f while IST data remains\n", ccd_time);
        if (ccd_time >= TIMELIMIT)
 		    fprintf(outgbias, "%15ld  %15.5f  %15.5f  %15.5f\n", b_est_cnt, b_est_sum/b_est_cnt*3600/PI*180.0);
 			error(msg);
        else
            error(msg) ;
        end
    end
    while (ccd_time - t_gyro) > 0.0001            
        pseudo_time = t_gyro ;
        if (pseudo_time < KF)  
            w = w0;
        end
        %rg = no_stars(pseudo_time, t_gyro, q0, w) ;
        rg = no_stars(q0, w) ;
        qp = q0;
        A = q_to_A(qp) ;		  
        M = zeros(3,3);
        for i=1:1:3 
            for j=1:1:3 
                for k=1:1:3
                    M(i,j) = M(i,j) + T_B(k,i) * A(k,j);
                end
            end
        end
        cnt0 = cnt0 + 1;
        if (cnt0 > 100.0)
            cctime = pseudo_time + Ctime ;
            consec = 0 ;
            QC1 = 0 ; atti_rms = -1.0 ;
        end
        if (cnt0 > 600.0)
            dt_merge = cnt0 ;
            ready_merge = 1 ;
            QC1 = 0 ; atti_rms = -2.0 ;
        end
        fprintf(out1b,"%15.6f  23 ", pseudo_time) ;  % No IST timetags */
        fprintf(out1b,"  %15.12f",qp) ;
        fprintf(out1b,"\n") ;
        ID_loca = '3' ; id_ = 23 ;
 		fprintf(qc,"%15.6f  %3d %6.2f\n", pseudo_time, QC1, atti_rms) ; stupid code
        fprintf(outn,"%15.6f  %2d  %2d  %2d  %2d  %8.3f\n", pseudo_time, id_, 0, 0, consec, 0.0) ;
		fprintf(outM,"%12.6f %3d\n%20.15f %20.15f %20.15f \n%20.15f %20.15f %20.15f \n% 20.15f %20.15f %20.15f \n", pseudo_time, QC1, M) ;
        [t_gyro, w, u, gread] = read_gyro(b_est_average, gyro_meas, gyro_count) ;
        
        if (gread == 1)
 		    fprintf(stdout,"No more gyro data. IST data may exist at %10.2f\n", t_gyro) ;
            if (ccd_time >= 86400.)
                fprintf(outgbias, "%15ld  %15.5f  %15.5f  %15.5f\n",  b_est_cnt, b_est_sum/b_est_cnt*3600/pi*180.0) ;
                error(msg) ;
            else
                error(msg) ;
            end
        end
    end % end of while (ccd_time == t_gyro) */
    %if (t_gyro < KF) for(i=0;i<3;i++) w[i] = w0[i]; 
    numlines = 0; % stupid code
    if ( (t_gyro - t_rms) >= DT_RMS) % rms calc. with 10 minute intervals */
        if (rmslnum > 0) 
            fprintf(outds," %12.5f %9d %12.8f %12.8f %12.8f\n", ccd_time, rmslnum,rmslsum/rmslnum*180.0/pi*3600) ;
		else
		    fprintf(outds," %12.5f %9d %12.8f %12.8f %12.8f\n",ccd_time, rmslnum, 0.0, 0.0, 0.0) ;
        end
        rmslsum = zeros(3,1) ;
		rmslnum = 0 ;
        
        if (rmszlnum > 0)
		    fprintf(outdq," %12.5f %9d %12.8f %12.8f %12.8f %12.8f\n", ccd_time, rmszlnum, rmszlsum/rmszlnum*180.0/3.1415*3600.0, rmszscalar/rmsznum ) ;
		else
		    fprintf(outdq," %12.5f %9d %12.8f %12.8f %12.8f %12.8f\n", ccd_time, rmszlnum, 0.0, 0.0, 0.0, 0.0) ;
		rmszlsum = zeros(3,1) ;
		rmszlscalar = 0.0 ;
		rmszlnum = 0 ;
		t_rms = t_rms + DT_RMS ;
        end
    end          
    if(ccd_time < KF) 
        w0 = w;
    end
	numlines = numlines + 1;
end
%end of while */
    
if (rmslnum > 0)
    fprintf(outds," %12.5f %9d %12.8f %12.8f %12.8f\n",ccd_time, rmslnum, rmslsum/rmslnum*180.0/3.1415*3600.0) ;
else  
    fprintf(outds," %12.5f %9d %12.8f %12.8f %12.8f\n",ccd_time, rmslnum, 0.0, 0.0, 0.0) ;
end
fprintf(outds," %12.5f %9d %12.8f %12.8f %12.8f\n",ccd_time, rmsnum, rmssum/rmsnum*180.0/3.1415*3600.0);
if (rmszlnum > 0)
    fprintf(outdq," %12.5f %9d %12.8f %12.8f %12.8f %12.8f\n",ccd_time, rmszlnum, rmszlsum/rmszlnum*180.0/3.1415*3600.0,rmszscalar/rmsznum ) ;
else  
    fprintf(outdq," %12.5f %9d %12.8f %12.8f %12.8f %12.8f\n",ccd_time, rmszlnum, 0.0, 0.0, 0.0, 0.0) ;
end
fprintf(outdq," %12.5f %9d %12.8f %12.8f %12.8f %12.8f\n",ccd_time, rmsznum,rmszsum/rmsznum*180.0/3.1415*3600.0,rmszscalar/rmsznum ) ;

fprintf(outgbias, "%15ld  %15.5f  %15.5f  %15.5f\n",b_est_cnt, b_est_sum/b_est_cnt*3600/3.1415*180.0) ;

fprintf(stdout, "End of PROGRAM\n") ; 

fprintf(stdout, "total_mstar = %15.0f  total_ided = %15.0f\n",total_mstar, total_ided) ;
fprintf(stdout, "c_pm = %8d  c_pm_ided = %8d  c_dm = %8d  c_dm_ided = %8d\n", c_pm, c_pm_ided, c_dm, c_dm_ided) ;

t1 = clock() ;
time(now2);
wallclock = difftime (now2, now) ;
fprintf(stderr, "It's now %s \n", ctime(now2)) ;
fprintf(stderr, "it took %.5f CPU secs.\n",(t1-t0)/CLOCKS_PER_SEC) ;
fprintf(stderr, "it took %.5f wall clock secs.\n", wallclock) ;
error(msg);