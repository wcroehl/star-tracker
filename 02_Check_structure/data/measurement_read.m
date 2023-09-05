clear all

% True initial attitue
q0_file = fopen('q_initial.txt', 'r'); 

t0 = fscanf(q0_file,'%g', 1);
q0 = fscanf(q0_file, '%g %g %g %g', 4)' ;

% Sensor measurement: Gyroscope
g_meas_file = fopen('g_rate_lrs.dat', 'r'); 

gyro_meas = [];
i = 1;
while( ~feof(g_meas_file) ) 
    gyro_time = fscanf(g_meas_file,'%g', 1);
    gyro_read = fscanf(g_meas_file, '%g %g %g', 3) ;
    gyro_meas = [gyro_meas; [gyro_time gyro_read']];
    i = i + 1
end

fclose(g_meas_file);

% INPUT : FOV read (Star camera) 
q_true_file = fopen('true_atti_lrs.dat', 'r'); 
q_true = [];
i = 1;
while( ~feof(q_true_file) ) 
    q_true_time = fscanf(q_true_file,'%g', 1);
    q_true_read = fscanf(q_true_file, '%g %g %g %g', 4) ;
    q_true = [q_true; [q_true_time q_true_read]];
    i = i + 1
end
fclose(q_true_file);

% INPUT : FOV read (Star camera) 
FOV_meas_file = fopen('FOVs.dat', 'r'); 
FOV_read = zeros(1,91);
FOVs_meas = [];
i = 1;
while( ~feof(FOV_meas_file) ) 
    FOV_time = fscanf(g_meas_file,'%g', 1);
    for j=1:1:91
        FOV_read(j) = fscanf(g_meas_file, '%g', 1) ;
    end
    FOVs_meas = [FOVs_meas; [FOV_time FOV_read]];
    i = i + 1
end

fclose(FOV_meas_file);

clearvars -except ...
t0 q0 ...
gyro_meas FOVs_meas

save('measurement_data.mat')
