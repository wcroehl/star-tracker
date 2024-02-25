% Dae Young Lee
% 8/22/16
%
% This is the wrapper file for testing the ICESat QUEST and running
% simulations to determine the overall attitude determination accuracy from
% star measurement vector.
%
% REVISION LOG
%
%   8/22/16 v1
%       initial version
%
% NOTES/INSTRUCTIONS
%

clear all
clc
close all

%% load the data (perfect measurements)
fid =fopen('QUEST.dat', 'w');
% directory where the data (.mat file) is:
dir_of_data = '.';
% filename. MODIFY THIS AS DESIRED:
flname = 'measurement.mat';

%% load the data
load([dir_of_data,'\',flname]);
% loaded variables are (Example):
%   vecTime          Times, seconds, m x 1, seconds elapsed

disp('Measurements loaded.')
q_ECI = q_true;

disp('Sim data loaded.');

%% Defne parameters to create measurments
ERR_POS = 0.3; % median error for Hipa. is 0.77 mas-> set to 0.1 
st_std = (6.5+ERR_POS)/3600.*pi/180.0;

%% Initial conditions and tuning parameters
Meas2Use_true = 1:1:length(q_true_time);
Meas2Use_vec = 1:1:length(vecTime);

c = 1; % to keep track of which index of gyro measurements we're on
disp('Running QUEST...');
for i = 1:length(Meas2Use_true)
    while vecTime(Meas2Use_vec(c)) == q_true_time(Meas2Use_true(i))
        % put the measured vectors togethter to pass to update function
        measVecs = vecMeas{c}';
        % build the corresponding reference vectors
        refVecs = vecRef{c}';
        refMags = magRef{c};

        [~,numVector] = size(refVecs);
        allSigs = repmat(st_std, numVector, 1);
        if (numVector >=3)
            [qhat_new, ~, dist_index, IFAIL] = q_method1(measVecs', refVecs', allSigs, refMags);
            fprintf(fid, '%15.9f %15.9f %15.9f %15.9f %15.9f \n', ...
                    vecTime(Meas2Use_vec(c)), ...
                    qhat_new(1), qhat_new(2), qhat_new(3), qhat_new(4));
        end
        c = c + 1;    % increment to the next measurement
        if (c <= length(Meas2Use_vec))
        else
            break;
        end
    end
    progressbar(i/(length(Meas2Use_true)))
end
fclose(fid);
quat_true_quest
