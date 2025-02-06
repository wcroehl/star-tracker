function ambi = sym_test(dist_tol, NEWBLI)

A = 1;
M = 2;
N = 3;

cata_ang = zeros(1,3);
%calculate possible cos angle
cata_ang(1) = sum(NEWBLI(A).L.*NEWBLI(M).L);
cata_ang(2) = sum(NEWBLI(A).L.*NEWBLI(N).L);
cata_ang(3) = sum(NEWBLI(M).L.*NEWBLI(N).L);

%change error value if exist
for i = 1:3
    if cata_ang(i) > 1
        cata_ang(i) = 1;
    elseif cata_ang(i) < -1
        cata_ang(i) = -1;
    end
    cata_ang(i) = acos(cata_ang(i));
end

%check whether cos angle diffrence less than tol(check whether stars do not create symetirc triangles)
if abs(cata_ang(1) - cata_ang(2)) < dist_tol || abs(cata_ang(1) - cata_ang(3)) < dist_tol || abs(cata_ang(2) - cata_ang(3)) < dist_tol
    ambi = 3; %symetric case
else
    ambi = 0; %non-symetric case
end