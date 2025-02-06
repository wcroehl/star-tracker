function [x, y, xmag, m_star] = obs_ang(nstar, x, y, xmag, m_star)

%save the first line value
x0 = m_star(1).L ; 
mag0 = m_star(1).mag ;
xx = x(1); 
yy = y(1);
mg = xmag(1);

%move up to save row
for j=2:1:nstar
    m_star(j-1).L   = m_star(j).L ; 
    m_star(j-1).mag = m_star(j).mag ;
    x(j-1) = x(j); 
    y(j-1) = y(j); 
    xmag(j-1) = xmag(j); 
end
 
%save the first line value to end
m_star(nstar).L   = x0; 
m_star(nstar).mag = mag0 ;
x(nstar) = xx ; 
y(nstar) = yy ; 
xmag(nstar) = mg ; 

