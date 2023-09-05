function supX = getSuperCross(v)
% Function takes the vector v and returns the skew symmetric matrix 
% representation of its cross product, i.e., the "super cross" matrix.
% Source: Aero 540, Aero 575, equation 3.149 in Crassidis, and
% http://en.wikipedia.org/wiki/Cross_product#Conversion_to_matrix_multiplication

supX = [0    -v(3)  v(2);
        v(3)  0    -v(1);
       -v(2)  v(1)  0];
   
end