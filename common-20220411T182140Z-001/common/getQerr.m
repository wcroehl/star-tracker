% getQerr returns the error quaternion for 2 quaternions or 2 arrays
% of quaternions.
% 
%   INPUTS:
%       q_est: 1 x 4 (for array modularity) estimated quaternion
%       q_true: 1 x 4 (for array modularity) true (or reference, generally)
%               quaternion
%
%   OUTPUTS:
%       q_err: Error between the two input quaternions
%
%   NOTE(s):
%       If the order of q_est and q_true is changed, the sign of the output
%       will flip.
%
%   SOURCE(s): 
%       - "Optimal Estimation of Dynamc Systems," Crassidis, Chapter 3
%
%% Author(s):
%   Alex Fox
%
%% REVISION LOG
%
%   8/29/12 -- v1
%       Created to ease plotting, it is also a useful function to have in
%       general
%
%% INSTRUCTIONS/NOTES
%
%

%%
function [q_err] = getQerr(q_est, q_true)

% pre-allocate the error quaternion array
q_err = zeros(length(q_est),4);

for i = 1:length(q_est)
    % Eq. 3.169 for inverse quaternion
    q_true_inv = [-q_true(i,1:3)'; q_true(i,4)];
    % Eq. 3.168 in source where q' = q_est, q = q_true^(-1)
    q_err(i,:) = [getXiMatrix(q_true_inv), q_true_inv]*q_est(i,:)';
end

end