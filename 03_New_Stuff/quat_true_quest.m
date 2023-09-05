%function quat_true_quest()

clear all
close all

%qtrue = fopen('true_atti_lrs.dat') ;
% qesti = fopen('ckg_1a.dat');
% qesti = fopen('QUEST.dat');

qesti = fopen('true_atti_lrs.dat','r') ;
qtrue = fopen('ckg_1a.dat','r');
delta = fopen('del_quest.dat','w') ;
bigD  = fopen('big_qs_delta.dat','w') ;

qtrue = fopen('QUEST.dat','r');

cat = fscanf(qtrue, '%s', 1);
ca  = fscanf(qtrue, '%s %s %s %s', 4);
cat = fscanf(qtrue, '%s', 1);
cqt = fscanf(qesti,'%s', 1)

delta_data = zeros(4,1);
bigD_data  = zeros(4,1);

while( ~feof(qesti) ) 
    cstk = fscanf(qesti,'%s', 1) ;
    q_true = fscanf(qesti, '%g %g %g %g', 4) ;
    t_q = str2num(cqt);
    
    if( abs(t_q - str2num(cat)) < 0.001 ) 
        q_esti = fscanf(qtrue, '%g %g %g %g', 4) ;
        t_a = str2num(cat);
        opp_sign = 0 ;
        for i=1:1:4
            if (q_esti(i)*q_true(i) < 0) opp_sign = opp_sign + 1; end
        end
        if (opp_sign >= 3) q_esti = -q_esti; end
        q_true = -q_true; 
        delqa11 = qcomp(q_true, q_esti)';
        delqa = getQerr(q_esti', q_true');
        delthea = 2 * delqa/norm(delqa);
        cat = fscanf(qtrue, '%s', 1);
        fprintf(delta, '%10.2f  %15.6f  %15.6f  %15.6f \n',...
                 t_q, delthea(1)*180./pi*3600,... 
                      delthea(2)*180./pi*3600, ...
                      delthea(3)*180./pi*3600) ;
        delta_data = [delta_data [t_q
                        delthea(1)*180./pi*3600
                        delthea(2)*180./pi*3600
                        delthea(3)*180./pi*3600] ];
                    
        if (abs(delthea(1)*180./pi*3600. > 20.0  || abs(delthea(2))*180./pi*3600. > 20.0))
           fprintf(bigD, '%10.2f  %15.6f  %15.6f  %15.6f \n', ...
                     t_q, delthea(1)*180./pi*3600, ...
	    	     delthea(2)*180./pi*3600, ...
		     delthea(3)*180./pi*3600) ;
            bigD_data = [bigD_data [t_q
                          delthea(1)*180./pi*3600
                          delthea(2)*180./pi*3600
                          delthea(3)*180./pi*3600] ];
        %    cat = fscanf(qtrue,'%s', 1) ;
        end
    end
    cqt = fscanf(qesti, '%s', 1)
end
fclose(qtrue) ;
fclose(qesti) ;
fclose(delta) ;
fclose(bigD) ;

figure;
subplot(3,1,1)
plot(delta_data(1,:),delta_data(2,:),'.');
axis([0 t_q -20 20]); grid on
xlabel('Time(sec)');ylabel('arcsec')
subplot(3,1,2)
plot(delta_data(1,:),delta_data(3,:),'.');
axis([0 t_q -20 20]); grid on
xlabel('Time(sec)');ylabel('arcsec')
subplot(3,1,3)
plot(delta_data(1,:),delta_data(4,:),'.');
axis([0 t_q -300 300]); grid on
xlabel('Time(sec)');ylabel('arcsec')
bigD(:,2:end)

%end
