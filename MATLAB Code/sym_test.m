function ambi = sym_test(dist_tol, NEWBLI)

A = 1;
M = 2;
N = 3;

cata_ang = zeros(1,3);

cata_ang(1) = sum(NEWBLI(A).L.*NEWBLI(M).L);
cata_ang(2) = sum(NEWBLI(A).L.*NEWBLI(N).L);
cata_ang(3) = sum(NEWBLI(M).L.*NEWBLI(N).L);

for i = 1:3
    if cata_ang(i) > 1
        cata_ang(i) = 1;
    elseif cata_ang(i) < -1
        cata_ang(i) = -1;
    end
    cata_ang(i) = acos(cata_ang(i));
end

if abs(cata_ang(1) - cata_ang(2)) < dist_tol || abs(cata_ang(1) - cata_ang(3)) < dist_tol || abs(cata_ang(2) - cata_ang(3)) < dist_tol
    ambi = 3;
else
    ambi = 0;
end