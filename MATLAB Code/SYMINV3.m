function [a, ifail] = SYMINV3(a)

ifail = 0 ;  %  construct truth table 
r = ones(3,1); % begin inversion 
p = ones(3,1);
q = [1/a(1,1) 1/a(2,2) 1/a(3,3)];

%  search for pivot 
for i=1:1:3 
    big = 0.0;
    for j=1:1:3
        test = abs(a(j,j)) ;
        if( (test-big) > 0 ) 
            if (r(j) < 0 ) 
                ifail = 1  ; 
                return ;
            end
            if (r(j) > 0 )
                big = test ; 
                k=j ; 
            end
        end
    end
     
    % preparation for elimination step 
	r(k)=0.0 ;
	a(k,k)=0.0 ;
    k1 = k+1 ;
    kp1= k+1 ;
    km1= k1-1 ;
    if (km1 < 1)
        ifail = 1 ;
        return ;
    end
	if (km1 > 1)
	     for j=1:1:km1 
             p(j) = a(j,k);
             q(j) = a(j,k) * q(k);
             if (r(j) < 0)
                 ifail = 1 ;
                 return ;
             end
             if (r(j) > 0)
                 q(j) = -q(j) ;
             end
             a(j,k) = 0.0 ;
         end
    end
    if ( (k1-4) > 0 )
        ifail = 1 ;
        return ;
    end
    if ( (k1-4) < 0 )
	     for j=kp1:1:3
             p(j) = a(k,j);
             if (r(j) < 0)
                 ifail = 1 ;
                 return;
             end
             if (r(j) == 0) 
                 p(j) = -p(j) ; 
             end
             q(j) = -a(k,j)*q(k);
             a(k,j) = 0.0 ;
         end
    end
    
    %  elimination proper 
    for j=1:1:3
        for k=j:1:3
          a(j,k) = a(j,k) + p(j)*q(k);
        end
    end
      
    %  replace lower half   
    for j=1:1:3
        for k=j:1:3
            a(k,j) = a(j,k) ;
        end
    end
end


% for(i=0;i<3;i++)  {        /*  search for pivot */
% 
% 	 k1 = k+1 ;
% 	 kp1=k+1 ;
% 	 km1=k1-1 ;
% 	 if (km1 < 0)  { *ifail = 1 ; return ; }
% 	 if (km1 > 0)
% 	     for (j=0 ; j<km1 ; j++) {
% 		 p[j] = a[j][k] ;
% 		 q[j] = a[j][k] * q[k] ;
% 		 if (r[j] < 0)  { *ifail = 1 ; return ; }
% 		 if (r[j] > 0)  q[j] = -q[j] ;
% 		 a[j][k] = 0. ;
% 	     }
% 	 if ( (k1-3) > 0 ) { *ifail = 1 ; return ; }
% 	 if ( (k1-3) < 0 )
% 	     for (j=kp1 ; j < 3 ; j++) {
% 		 p[j] = a[k][j] ;
% 		 if (r[j] < 0) { *ifail = 1 ; return ; }
% 		 if (r[j] == 0)  p[j] = -p[j] ;
% 		 q[j] = -a[k][j]*q[k] ;
% 		 a[k][j] = 0.0 ;
% 	 }
% }

