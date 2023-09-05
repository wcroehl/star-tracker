function [rg, q] = no_stars(q, w)  

dq=zeros(4,1);% added code, 7/14/23
q = q/norm(q);
            
dq(1) =( w(3)*q(2)-w(2)*q(3)+w(1)*q(4))/2.0/10.0 ; 
dq(2) =(-w(3)*q(1)+w(1)*q(3)+w(2)*q(4))/2.0/10.0 ;
dq(3) =( w(2)*q(1)-w(1)*q(2)+w(3)*q(4))/2.0/10.0 ;
dq(4) =(-w(1)*q(1)-w(2)*q(2)-w(3)*q(3))/2.0/10.0 ;

q = q + dq; 
q = q/norm(q);
rg = 2;

end