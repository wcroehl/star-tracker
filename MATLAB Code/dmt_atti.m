% void dmt_atti(w, q0, Mp)  /* propagate attitude using the gyro meas. */ 
function Mp = dmt_atti(w, q0,T_B)  % propagate attitude using the gyro meas. % 
% double w[3], q0[4], Mp[3][3] ;
%  {
%  int     i, j, k ;
%  double  w_mag, n_hat[3], cwt, swt, qp[4], Ap[3][3] ;

T_STEP_GYRO = 0.1;
%T_B = [1 0 0;0 1 0;0 0 1]; remove, 09/23/23

 w_mag = sqrt(w(1)*w(1)+w(2)*w(2)+w(3)*w(3)) ;
 for(i=1:3)  n_hat(i) = w(i)/w_mag ; end

 cwt = cos(w_mag*T_STEP_GYRO/2.0) ;  % Originally multiplied by 10 %
 swt = sin(w_mag*T_STEP_GYRO/2.0) ;  % Originally multiplied by 10 %

 qp(1) = cwt*q0(1) + swt*(n_hat(3)*q0(2) - n_hat(2)*q0(3) + n_hat(1)*q0(4)) ;
 qp(2) = cwt*q0(2) + swt*(-n_hat(3)*q0(1) + n_hat(1)*q0(3) + n_hat(2)*q0(4)) ;
 qp(3) = cwt*q0(3) + swt*(n_hat(2)*q0(1) - n_hat(1)*q0(2) + n_hat(3)*q0(4)) ;
 qp(4) = cwt*q0(4) + swt*(-n_hat(1)*q0(1) - n_hat(2)*q0(2) - n_hat(3)*q0(3))  ;

 Ap = q_to_A(qp) ;

 Mp = zeros(3,3) ;
 for(i=1:3) 
     for(j=1:3) 
         for(k=1:3)
             Mp(i,j) = Mp(i,j) + T_B(k,i) * Ap(k,j) ; % Mp = Ap for IST %
         end
     end
 end
end