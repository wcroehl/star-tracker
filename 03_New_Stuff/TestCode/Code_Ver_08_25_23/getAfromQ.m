
% function to calculate a rotation matrix A from quaternions
%
% Since the quaternions represent attitude of the body frame relative to
% the inertial frame, the resulting attitude matrix A is equuivalent to the
% Orientation matrix O_BI, that is, the orientation of the body frame
% relative to the inertial frame
%
% Author(s)
%   John Springmann
%   8/2/12
%
% The source of the equations is Section 3.7.1 of Crassidis' "Optimal
% Estimation of Dynamc Systems"

function A = getAfromQ(q)
% INPUT q is 4 x 1 quation with elements 1-3 being the vector component

    q1 = q(1);
    q2 = q(2);
    q3 = q(3);
    q4 = q(4);

    Xi = getXiMatrix(q);
    Phi = [q4*eye(3) - getSuperCross(q(1:3)); -[q1, q2, q3]];
    
    A = Xi' * Phi;
end