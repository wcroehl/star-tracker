clear,clc
load('measurement_data.mat');
cdata = fopen("C:\Users\Public\research\EKF_c\data\FOVs2.txt");
cdata2 = fscanf(cdata, "%lf");
for i = 0 : 5162428
    m = mod(i,91)+1;
    d = floor(i/91)+1;
    daduh(d,m) = cdata2(i+1);
end
for x = 1:91
    for y = 1:56728
        if(FOVs_meas(y,x) ~= daduh(y,x))
            x
            y
        end
    end
end