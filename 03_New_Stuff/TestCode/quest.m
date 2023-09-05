function Qest = quest(W, V, a)

N = length(a);
B = zeros(3,3);
for i=1:N,
    B = B + a(i)*(W(:,i)*V(:,i)');
end

S   = B+B';
rho = trace(B);

Z = zeros(3,1);
for i=1:N,
    Z = Z + a(i)*cross(W(:,i),V(:,i)); %this is the only line I wasn't sure about - for now
                                    %I am assuming this is correct
end

K = [ S-rho*eye(3) Z
      Z'           rho ];
  % in Forbes textbook, rho = k22, Z = k12

%[E,D] = eig(K);          % Find the eigenvector for largest eigenvalue of K
%[a,b] = max(diag(D));

%qq      = E(:,b);         % Compute the unit quaternion
%q      = qq/norm(qq);

%guess the largest eigenvalue - use sum of weights, then input into
%Newton's method to get a better largest eigenvalue guess
%using the assumption of equal weights above, the guessed e-value is simple
g = 1;
eguess = newton(rho, S, Z, g);

%QUEST method
%compute coefficients
alpha = eguess^2 - rho^2 + trace(adj(S)); %need classical adjoint of S
beta = eguess - rho;
gamma = (eguess+rho)*alpha - det(S);

%compute x
x = (alpha*eye(3) + beta*S + S^2)*Z;

%plug all into formula to get quaternion estimate q2
Qest = 1/(sqrt(gamma^2+x'*x))*[x;gamma];

end

function B = adj(A)
%computes the adjugate matrix (ie, the classical adjoint) of A
C = cof(A);
B = C';

end

function C = cof(A)
%computes the cofactor matrix of A
%A assumed to be square
[m,n] = size(A);
C = zeros(m,m);

for i=1:1:m
   for j=1:1:m
       M = A;
       C(i,j) = (-1)^(i+j)*det(M);
   end
end

end

function bestguess = newton(rho, S, Z, guess)
%coded for use with QUEST
%uses Newton-Raphson method to compute best guess for largest eigenvalue
%the tolerance eps can be changed as needed, though diff must be greater
%than eps for the first iteration
a = rho^2-trace(adj(S));
b = rho^2+Z'*Z;
c = det(S) + Z'*S*Z;
d = Z'*S^2*Z;
eps = 1e-8; %a tolerance
diff = 90; %arbitrarily large for first iteration
x1 = guess;
while (diff >= eps)
   x2 = x1 - f(a, b, c, d, rho, x1)/fprime(a, b, c, d, rho, x1);
   diff = abs(x2-x1);
   x1 = x2;
end
bestguess = x1;
end

function y = f(a, b, c, d, rho, x)
%polynomial equation for finding largest eigenvalue estimate; this is
%derived from the q-method eigenproblem
%a, b, c, d are computed in newton() and brought in here
y = x^4 - (a+b)*x^2 - c*x + (a*b + c*rho - d);
end

function y = fprime(a, b, c, d, rho, x)
%polynomial derivative with respect to x (the eigenvalue estimate)
y = 4*x^3 - 2*(a+b)*x - c;
end