%int  cuvf(t_new, t_old, t_elapse, qq, q, P, w, b, 
%          cnt, zsum, znum, ccdstep)
%           function [] = curvf(SIG_ARW, SIG_RRW )
function [zsum, znum] = cuvf(t_new, t_old, w)

end
clear all
SIG_ARW = 1.0e-8;
sigv = SIG_ARW;
SIG_RRW = 1.0e-8;
sigu = SIG_RRW;
ERR_POS = 0.3;          
s_a = SIG_ARW*SIG_ARW ;
s_r = SIG_RRW*SIG_RRW ;
std_dev = (6.5+ERR_POS)/3600.*pi/180.0 ;
R = diag( [std_dev*std_dev, std_dev*std_dev, std_dev*std_dev]) ;
  %for(i=0;i<3;i++) P[i][i] = 1.0e-6 ;
  %for(i=3;i<6;i++) P[i][i] = 1.0e-8 ;
  P = diag([1.0e-6 1.0e-6 1.0e-6 1.0e-8 1.0e-8 1.0e-8]);
  P0 = P;


w_here = [ 0.3 0.1 0.2]';
w = w_here;
Omega = diag([w(1), w(2), w(3)]);

 deltaT = 0.01;
del_t = deltaT;
ddel_t = del_t*del_t ;

% % eqn 7.55:
% psi = sin( 1/2* norm(w_here) * deltaT) * w_here / norm(w_here);
% % supercross
% psi_sup = getSuperCross(psi);
% 
% % eqn 7.54:
% OMEGA = [ cos( 1/2*norm(w_here)*deltaT ) * eye(3) - psi_sup, ...
%               psi; ...
%               -psi', ...
%               cos( 1/2*norm(w_here)*deltaT )]
%           
% % Propagate covariance: --------------------------------------------
% 
% % super-cross of estimated angular rate
% w_here_sup = getSuperCross( w_here);
% 
% % elements of the state transision matrix
% phi11 = eye(3) - w_here_sup* sin(norm(w_here)*deltaT)/norm(w_here) + w_here_sup^2 * ( 1 - cos( norm(w_here)*deltaT ))/norm(w_here)^2 ;
% phi12 = w_here_sup * (1 - cos(norm(w_here)*deltaT) )/norm(w_here)^2 - eye(3)*deltaT - w_here_sup^2*( norm(w_here)*deltaT - sin(norm(w_here)*deltaT))/norm(w_here)^3;
% phi21 = zeros(3,3);
% phi22 = eye(3);
% % build state transision matrix:
% phi1 = [phi11 phi12; phi21 phi22]
% 
% % process covariance matrix
% Q = [( sigv^2*deltaT + 1/3*sigu^2*deltaT^3 )*eye(3), ...
%     (1/2*sigu^2*deltaT^2)*eye(3); ...
%     (1/2*sigu^2*deltaT^2)*eye(3), ...
%     (sigu^2*deltaT)*eye(3)];

% Tune:
% Q = 20.*Q;
%Q = tuneFac.*Q; % new input

% propogate:
% G = [-eye(3) zeros(3,3); zeros(3,3) eye(3)];
% Pnext = phi*P0*phi' + G*Q*G';
% 
% q = qnext;
% P = Pnext;          
%           

        w_mag = norm(w) ; 
        n_hat = w/w_mag ; 
        phi = w_mag * deltaT;   
        
        % State transition matrix */
        n33 = antisym3(n_hat) ;
        n33sq = zeros(3,3);
        n33sq = n33*n33;
        for i=1:1:3
            for j=1:1:3
                s_n33(i,j) = n33(i,j)*sin(phi) ;
                c_n33sq(i,j) = n33sq(i,j) * (1-cos(phi)) ;
                transit(i,j) = s_n33(i,j) + c_n33sq(i,j) ;
            end
            transit(i,i) = transit(i,i)+ 1.0 ;
        end
        transit1(1:3,1:3) = eye(3) + sin(phi).*n33 + (1-cos(phi)).*n33sq;
      a_ = deltaT ;   ksai = ( 1 - cos(phi) )/ phi ;
      b_ = a_ * ksai ;          eta  = 1 - sin(phi)/phi ;
      c_ = a_ * eta ;

      for i=1:1:3
          for j=4:1:6
              transit(i,j) = n33(i,j-3)*b_ + n33sq(i,j-3)*c_ ;
          end
          transit(i,i+3) = transit(i,i+3) + a_ ;
      end
      
      transit1(1:3,4:6) = b_.*n33 + c_.*n33sq + a_*eye(3,3); 

      for i=4:1:6
          for j=1:1:6
              transit(i,j) = 0.0 ;
          end
          transit(i,i) = 1.0 ;
      end
            transit1(4:6,1:3) = zeros(3,3)
      transit1(4:6,4:6) = eye(3) ;

      % P propagated 
      P_propa = transit*P0*transit';
      
% 
%       /* Q_k matrix */
% 
%       for(i=0;i<6;i++) for(j=0;j<6;j++) Q_k[i][j] = 0.0 ;
% 
%       /* 1 to 3 rows */
%       for(i=0;i<3;i++)
%                 Q_k[i][i] = (s_a + s_r*ddel_t/3.0)*del_t ;
%       for(i=0;i<3;i++)
%                 Q_k[i][i+3] = -s_r*ddel_t/2.0 ;
% 
%       /* 1 to 3 columns */
%       for(i=0;i<3;i++)
%                 Q_k[i+3][i] = -s_r*ddel_t/2.0 ;
% 
%       /* diagonal terms */
%       for(i=0;i<3;i++)
%                 Q_k[i+3][i+3] = s_r*del_t ;  /* ddel_t ? */
% 
%       for(i=0;i<6;i++) for(j=0;j<6;j++) P_propa[i][j] += Q_k[i][j] ; 
%       
%       for(i=0;i<6;i++) for(j=0;j<6;j++)
%               P_pro[i][j] = (P_propa[i][j]+P_propa[j][i])/2.;
% %         
% %         % Q_k matrix
% %         Q_k = zeros(6,6);
% %         
% %         % 1 to 3 rows 
% %         Q_k_temp = (s_a + s_r*ddel_t/3.0)*del_t ;
% %         Q_k(1:3,1:3) = diag([Q_k_temp Q_k_temp Q_k_temp]);
% %         
% %         Q_k_temp = -s_r*ddel_t/2.0 ;
% %         Q_k(1:3,4:6) = diag([Q_k_temp Q_k_temp Q_k_temp]);
% %         
% %         % 1 to 3 columns 
% %         Q_k(4:6,1:3) = diag([Q_k_temp Q_k_temp Q_k_temp]);
% %         
% %         % diagonal terms 
% %         Q_k_temp = s_r*del_t;
% %         Q_k(4:6,4:6) = diag([Q_k_temp Q_k_temp Q_k_temp]);
% %         
% %         P_propa = Q_k;  
% %         for i=1:1:6
% %             for j=1:1:6
% %                 P_pro(i,j) = (P_propa(i,j)+P_propa(j,i))/2.;
% %             end
% %         end
% %         
% %         
% % 
% %      }  
% % 
% %   if(abs(t_new - t_old) < 0.0001 )
% %       P_pro = P; 
% %       q_propa = q ;
% %   end  
% %       
% %    
% % % -----  update  ------ 
% % A_pro = getAfromQ(q_propa);
% % W_pro = A_pro* crf(ic).vec;
% % W_anti = getSuperCross(W_pro);
% % H = zeros(3,6);
% % H(1:3,1:3) = W_anti;
% % PHT = zeros(6,3);
% % PHT = P_pro*H;
% % HPHT = H*PHT; 
% % 
% % % matrix B 
% % B = HPHT + R;
% %         
% % %B(1,2) = B(2,1) = (B(1,2)+B(2,1))/2.0 ;
% % %B(1,3) = B(3,1) = (B(1,3)+B(3,1))/2.0 ;
% % %B(2,3) = B(3,2) = (B(2,3)+B(3,2))/2.0 ;
% % 
% % % SYMINV3(B,&ifail) ;
% % % if(ifail == 1)
% % %     {
% % %    printf(" Inverse of B matirx is singular\n") ;
% % %    exit(1) ; }
% % % B[0][1] = B[1][0] = (B[0][1]+B[1][0])/2.0 ;
% % % B[0][2] = B[2][0] = (B[0][2]+B[2][0])/2.0 ;
% % % B[1][2] = B[2][1] = (B[1][2]+B[2][1])/2.0 ;
% % 
% % K = PHT*B;
% % I_KH = eye(6)-K*H; %(6,3)(3,6)
% % P = I_KH*P_pro;
% % P_new = P*I_KH;
% % KR = K*R;
% % P_pro = KR*K;
% % 
% % % for(i=0;i<6;i++) if(P_pro[i][i] < 0.0 ) P_pro[i][i] = 0.0 ;
% % % for(i=0;i<6;i++) for(j=0;j<6;j++)
% % %         P_new[i][j] += P_pro[i][j] ;
% % % for(i=0;i<6;i++) for(j=i;j<6;j++)
% % %         P[i][j] = P[j][i] = (P_new[i][j] + P_new[j][i])/2. ;
% % %         
% % % for(i=0;i<6;i++) if(P[i][i] < 0.0 ) P[i][i] = 0.0 ;
% % 
% % % dx update
% % %for(i=0;i<3;i++) z[i] = b_star[ic]->L[i] - W_pro[i] ;
% % z = b_star[ic].vec - W_pro;    
% % %for(i=0;i<3;i++) zsum[i] += sqrt(z[i]*z[i]) ;
% % %                 (*znum)++ ;
% % zsum = zsum + z.*z;                 
% % 
% % dx = zeros(6,1);
% % dx = dx + K*z % K = 6x3        
% % 
% % temp = dx(1:3)'*dx(1:3); 
% % coef = 1/sqrt(1.+ temp/4.) ;   
% %   
% % dq(1:3) = -coef *dx/2 ;  %/* +coef */ 
% % dq(4)   = coef ;
% %   
% % qcomp(dq,q_propa,q_new) ; 
% % 
% % q_new = q_new/norm(q_new) ;
% % b = b + dx(4:6);
% % q = q_new ;
% % t_old = t_new ;
% % 
% % }
% % 
% % return(0) ;
% % 
% % } 
% % 
% % 
