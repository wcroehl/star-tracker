function  [q, lambda] = fastq(S, Z, rho) 
% calculate optimal q 
% N-R Method to get the largest lambda with initial value  1 

a = rho^2 - trace(getAdjoint(S)) ; % compare with lines 100 to 115 of q_method.c, 7/17/23
b = rho^2+Z'*Z; %check whether rho is the smae as sig in C, 7/17/23
c = det(S) + Z'*S*Z; 
d = Z'*S^2*Z;
%equation from pg.473 of "Space Dynamics ..." book 

eps = 1e-8;    %a tolerance
diff = 1;      %arbitrarily large for first iteration
lambda = 1.0 ; % initial guess 
while (diff >= eps) % comapre the code around here with lines 119 and rtnewt function of C, 7/17/23
   x2 = lambda - func(a, b, c, d, rho, lambda)/fprime(a, b, c, d, rho, lambda); 
   diff = abs(x2-lambda);
   lambda = x2;
end

% QUEST method
% Calculate eigenvector corres. to lambda */
alpha = lambda^2 - rho^2 + trace(getAdjoint(S)) ;
beta  = lambda - rho ;
gamma = (lambda+rho)*alpha - det(S) ; 
%equations from pg.472 of "Space Dynamics ..." book 

%compute x
x = (alpha*eye(3) + beta*S + S^2)*Z; % comapre the calcualtions with lines 128 to 140 of q_method.c, 7/17/23
                                     %equation from pg.472 of "Space Dynamics ..." book 

%plug all into formula to get quaternion estimate q2

q = 1/(sqrt(gamma^2+x'*x))*[x;gamma]; %equation from eq 25.17 in pg.472 of "Space Dynamics..." book

end

function y = func(a, b, c, d, rho, x)
%polynomial equation for finding largest eigenvalue estimate; this is
%derived from the q-method eigenproblem
%a, b, c, d are computed in newton() and brought in here
y = x^4 - (a+b)*x^2 - c*x + (a*b + c*rho - d); %equation from eq 25.20 in pg.473 of "Space Dynamics ..." book 
end

function dy = fprime(a, b, c, d, rho, x)
%polynomial derivative with respect to x (the eigenvalue estimate)
dy = 4*x^3 - 2*(a+b)*x - c; % dertive equation from eq 25.20 in pg.473 of "Space Dynamics ..." book 
end

