clear all, clear clc
%Q_method  test code 

%load and make necessary files
outPTk2 = fopen('outpt2k.dat', 'w');
outPTk3 = fopen('outpt3k.dat', 'w');
out1a = fopen('ckg_1ak.dat', 'w');
out2 = fopen('ckg_2k.dat', 'w');
n_test=load('On_QTest.dat');
t_test=load('Ot_QTest.dat');
V_test=load('OV_QTest.dat');
W_test=load('OW_QTest.dat');
mag_test=load('Omag_QTest.dat');

%set the data for q_method test
l=length(n_test);
S=1;
while S <= l
    nstar=n_test(S);
    ccd_time=t_test(S);
    d_i=0;
    sign1=0;
    nsig=0;
    [rowV, colV] = find(V_test == ccd_time);
    [rowW, colW] = find(W_test == ccd_time);
    if  nstar == 4
       V=V_test(rowV:rowV+3,1:3);
       V(rowV:rowV+3,4)=mag_test(S,1:4);
       W=W_test(rowW:rowW+3,1:4);
    elseif nstar == 5
       V=V_test(rowV:rowV+4,1:3);
       V(rowV:rowV+4,4)=mag_test(S,1:5);
       W=W_test(rowW:rowW+4,1:4);
    else 
       V=V_test(rowV:rowV+5,1:3);
       V(rowV:rowV+5,4)=mag_test(S,1:6);
       W=W_test(rowW:rowW+5,1:4);
    end
    
    %q_method test
    [qq, P_theta, d_i, rtn] = q_method(W, V, nstar, d_i, outPTk2, outPTk3, ccd_time);
    
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
   S=S+1;
 end