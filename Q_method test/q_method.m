function  [q_est, P, dist_index, IFAIL] = q_method(W, V, nstar, dist_index, outPT2, outPT3, t) 
% add outPt's and t
ERR_POS = 0.3;

[nstar, ~] = size(W);
dist_index = 0;
wgt = zeros(nstar,1);
for i=1:1:nstar-1
    dist_index = dist_index + W(i,1)*W(i+1,1) + W(i,2)*W(i+1,2) + W(i,3)*W(i+1,3);
    %wgt(i) = 1/nstar;
end

%suggestion code
for i=1:nstar
    wgt(i) = 1/nstar;
end

dist_index = dist_index + W(nstar,1)*W(1,1) + W(nstar,2)*W(1,2) + W(nstar,3)*W(1,3);

%B = W'.*diag(wgt)* V;
B = zeros(3);
for i = 1:nstar
    for j = 1:3
        for k = 1:3
        B(j,k) = B(j,k) + wgt(i) * W(i,j) * V(i,k);
        end
    end
end
rho = trace(B); 
S = B + B';
Z = zeros(3,1);

for i=1:1:nstar
    %sp = getSuperCross(W(i,:))
    %Z = Z + wgt(i)*sp.*V(i,1:3)';
    Z(1) = Z(1) + wgt(i)*(W(i,2)*V(i,3)- W(i,3)*V(i,2));
    Z(2) = Z(2) + wgt(i)*(W(i,3)*V(i,1)- W(i,1)*V(i,3));
    Z(3) = Z(3) + wgt(i)*(W(i,1)*V(i,2)- W(i,2)*V(i,1));
end
lambda = 0;   
[q_est, lambda] = fastq(S, Z, rho);

A = q_to_A(q_est);  
v = V(1:3,1:3);
d = diag(wgt);
D = d(1:3,1:3);
P = zeros(3);

for i = 1:nstar
    for j = 1:3
        AV(j) = A(j,1)*V(i,1) + A(j,2)*V(i,2) + A(j,3)*V(i,3);
    end
    for j = 1:3
        for k = 1:3
            P(j,k) = P(j,k) + wgt(i)*AV(j)*AV(k);
        end
    end
end

for i = 1:3
    P(i,i) = 1 - P(i,i);
end

%add fprintf, 07/28/23
fprintf(outPT2,'%15.6f \n%.15f %.15f %.15f \n%.15f %.15f %.15f \n%.15f %.15f %.15f \n',t, P);

dev_tot = 0.0 ;
mag = V(:,4);
for i=1:1:nstar 
   if (mag(i) >= 2.0) 
       std_dev = (4.5+0.7*(mag(i)-2.)+ERR_POS)/3600.*pi/180.;
   end
   if (mag(i) <  2.0) 
       std_dev = (4.5+ERR_POS)/3600.*pi/180. ; 
   end
   dev_tot = dev_tot + 1/(std_dev*std_dev) ;
end

% P(2,1) = (P(1,2)+P(2,1))/2.0 ;
% P(1,2) = P(2,1);
% P(3,1) = (P(1,3)+P(3,1))/2.0 ;
% P(1,3) = P(3,1);
% P(3,2) = (P(2,3)+P(3,2))/2.0 ;
% P(2,3) = P(3,2);
% P 
P = (P + P')/2;
P = inv(P);

if det(P) > 1e-6
    IFAIL = 0;
else
    IFAIL = 1;
end

%[P, IFAIL] = SYMINV3(P);


% if (IFAIL == 1)   
%     error('Matrix inversion error at q_method()') ;
%     return;
% IFAIL = 0;

% P(2,1) = (P(1,2)+P(2,1))/2.0 ;
% P(3,1) = (P(1,3)+P(3,1))/2.0 ;
% P(3,2) = (P(2,3)+P(3,2))/2.0 ;
% P(1,3) = P(3,1);
% P(1,2) = P(2,1);
% P(2,3) = P(3,2);
P = (P + P')/2;
%     for i=1:1:3
%         for j=1:1:3
%             P(i,j) = P(i,j)/dev_tot ; % P_dthe_dthe = 4*P_dQ_dQ */
%         end
%     end
%     
    P = P/dev_tot;
%add fprintf, 07/28/23
fprintf(outPT3, '%15.6f \n%.15f %.15f %.15f \n%.15f %.15f %.15f \n%.15f %.15f %.15f \n',t, P);

P_range = zeros(3,1);    
P_range(1) = sqrt(4*P(1,1))*180/pi*3600;
P_range(2) = sqrt(4*P(2,2))*180/pi*3600;
P_range(3) = sqrt(4*P(3,3))*180/pi*3600;    
        
end 