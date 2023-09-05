function F = getXiMatrix(q)
% John Springmann
% 8/2/12
%
% equation 3.155a in Crassidis
%
% INPUT q is 4 x 1 quation with elements 1-3 being the vector component

q1 = q(1);
q2 = q(2);
q3 = q(3);
q4 = q(4);

super = getSuperCross(q(1:3));

F = [q4*eye(3) + super; -[q1, q2 q3]];

end