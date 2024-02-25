function quat_true_quest()

qtrue = fopen('true_atti_lrs.dat') ;
qesti = fopen('ckg_1a.dat');
%qesti = fopen('QUEST.dat');

cat = fscanf(qesti, '%s', 1);
ca  = fscanf(qesti, '%s %s %s %s', 4);

cat = fscanf(qesti, '%s', 1);
cqt = fscanf(qtrue, '%s', 1);
delta = zeros(4,1);
bigD  = zeros(4,1);

while( ~feof(qtrue) ) 
    cstk = fscanf(qtrue,'%s', 1) ;
    q_true = fscanf(qtrue, '%g %g %g %g', 4) ;
    t_q = str2num(cqt);
    if( abs(t_q - str2num(cat)) < 0.001 ) 
        q_esti = fscanf(qesti, '%g %g %g %g', 4) ;
        t_a = str2num(cat);
        opp_sign = 0 ;
        for i=1:1:4
            if (q_esti(i)*q_true(i) < 0) opp_sign = opp_sign + 1; end
        end
        if (opp_sign >= 3) q_esti = -q_esti; end
        delqa = getQerr(q_esti', q_true');
        delthea = 2 * delqa/norm(delqa);
        cat = fscanf(qesti, '%s', 1);
        delta = [delta [t_q
                        delthea(1)*180./pi*3600
                        delthea(2)*180./pi*3600
                        delthea(3)*180./pi*3600] ];
                    
        if (abs(delthea(1)*180./pi*3600. > 20.0  || abs(delthea(2))*180./pi*3600. > 20.0))
            bigD = [bigD [t_q
                          delthea(1)*180./pi*3600
                          delthea(2)*180./pi*3600
                          delthea(3)*180./pi*3600] ];
        end
    end
    cqt = fscanf(qtrue, '%s', 1)
end
fclose(qtrue);
fclose(qesti);

figure;
subplot(3,1,1)
plot(delta(1,:),delta(2,:),'.');
axis([0 t_q -20 20]); grid on
xlabel('Time(sec)');ylabel('arcsec')
subplot(3,1,2)
plot(delta(1,:),delta(3,:),'.');
axis([0 t_q -20 20]); grid on
xlabel('Time(sec)');ylabel('arcsec')
subplot(3,1,3)
plot(delta(1,:),delta(4,:),'.');
axis([0 t_q -300 300]); grid on
xlabel('Time(sec)');ylabel('arcsec')
bigD(:,2:end)

end
