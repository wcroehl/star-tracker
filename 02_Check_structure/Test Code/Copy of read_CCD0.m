function [t, x, y, xmag, m_star, nstar] = read_CCD0(FOVs_meas, Ns, Counter)  
%Charles

t = FOVs_meas(Counter,1);

nstar = 0 ;
for i=1:1:Ns
    cx(i)  = FOVs_meas(Counter,2+3*(i-1)); 
    cy(i)  = FOVs_meas(Counter,3+3*(i-1)); 
    mag(i) = FOVs_meas(Counter,4+3*(i-1)); 
    m_star(i).L = zeros(3,1);
    m_star(i).mag = mag(i);
    phi(i) = cx(i);
    lambda(i) = cy(i);
    if (mag(i) >= -2.0 && mag(i) < 7.5)  
        nstar = nstar + 1;
    end
end

for i=1:1:nstar
    x(i) = phi(i);
    y(i) = lambda(i);
    xmag(i) = mag(i);
    phi(i) = cx(i)*pi/180.0 ; % degree to radian
    lambda(i) = cy(i)*pi/180.0 ; % degree to radian 
    s(1) = phi(i);    % tan(phi) in original code */
    s(2) = lambda(i); % tan(lambda) in original code */
    m_star(i).L(3) = 1./sqrt(1+s(1)*s(1)+s(2)*s(2));
    m_star(i).L(1) = s(1)*m_star(i).L(3);
    m_star(i).L(2) = s(2)*m_star(i).L(3);
    m_star(i).mag = mag(i) ;
end
 
for i=1:1:nstar
    for j=i+1:1:nstar % bright star first 
        if(m_star(j).mag < m_star(i).mag)
            temporary.mag = m_star(i).mag;
            temporary.L   = m_star(i).L;
            tmp_x = x(i) ;
            tmp_y = y(i) ;
            tmp_mag = xmag(i); 
            m_star(i).mag = m_star(j).mag;
            m_star(i).L = m_star(j).L;
            x(i) = x(j) ; y(i) = y(j) ; xmag(i) = xmag(j) ; 
            m_star(j).mag = temporary.mag;
            m_star(j).L = temporary.L;
            x(j) = tmp_x ; y(j) = tmp_y ; xmag(j) = tmp_mag ; 
        end
    end
end

    