function R = getRotation(axis, angleRad)
% Function produces rotation matrix R for a rotation of angleRad radians
% about the axis specified.  Source: Aero 540 notes which match the
% equations cited on
% http://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula#Conversion_to_rotation_matrix

% Convert axis to column vector if it is a row vector
if sum(size(axis) == [1 3]) == 2
    axis = axis';
end


% Make rotation matrix for given axis and angle
R = eye(3) + getSuperCross(axis)*sin(angleRad) + (1 - cos(angleRad))*(axis*axis' - eye(3));

end