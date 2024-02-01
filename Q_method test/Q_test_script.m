clear all, clear clc
%Q_method  test code 

%load and make necessary files
outPTk2 = fopen('outptk2.dat', 'w');
outPTk3 = fopen('outptk3.dat', 'w');
n_test=load('On_QTest.dat');
t_test=load('Ot_QTest.dat');
V_test=load('OV_QTest.dat');
W_test=load('OW_QTest.dat');

%set the data for q_method test
l=length(n_test);
i=1;
while i <= l
    nstar=n_test(i);
    ccd_time=t_test(i);
    d_i=0;
    [rowV, colV] = find(V_test == ccd_time);
    [rowW, colW] = find(W_test == ccd_time);
    if  nstar == 4
       V=V_test(rowV:rowV+3,1:4);
       W=W_test(rowW:rowW+3,1:4);
    elseif nstar == 5
       V=V_test(rowV:rowV+4,1:4);
       W=W_test(rowW:rowW+4,1:4);
    else 
       V=V_test(rowV:rowV+5,1:4);
       W=W_test(rowW:rowW+5,1:4);
    end
    
    %q_method test
    [qq, P_theta, d_i, rtn] = q_method(W, V, nstar, d_i, outPTk2, outPTk3, ccd_time);
    i=i+1;
 end