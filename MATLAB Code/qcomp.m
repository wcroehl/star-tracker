function q_new = qcomp(a,b)
q_new = zeros(4,1);
q_new(1) = a(4)*b(1) + a(3)*b(2) - a(2)*b(3) + a(1)*b(4);
q_new(2) = -a(3)*b(1) + a(4)*b(2) + a(1)*b(3) + a(2)*b(4);
q_new(3) = a(2)*b(1) - a(1)*b(2) + a(4)*b(3) + a(3)*b(4);
q_new(4) = -a(1)*b(1) - a(2)*b(2) - a(3)*b(3) + a(4)*b(4);
end