function [t_gyro, w, u, end_indicator, counter] = read_gyro(b, gyro_meas, counter)
w=zeros(3,1);
if length(gyro_meas) > counter % This is for the end of the data.
    t_gyro = gyro_meas(counter,1); % time
    u_arcsec = gyro_meas(counter,2:4); % x, y, z gyro
    u = u_arcsec/3600*pi/180; % Arcsec to radian..
    w = u - b; % Gyro has bias but bias correction happens here.
    end_indicator = 0;
else
    t_gyro = [];
    w = [];
    u = [];
    end_indicator = 1;
end
counter = counter + 1;
%% Original C codes

% int read_gyro(t_gyro, u, b, w, gyro)
% FILE      *gyro ;
% double    *t_gyro, u[3], b[3], w[3] ;
%   {
%   char cgt[16], cx[16], cy[16], cz[16] ; 
%   int  i ;
% 
%   if (fscanf(gyro, "%s", cgt) != EOF)
%     {
%     *t_gyro = atof(cgt) ;
%     fscanf(gyro, "%s %s %s", cx, cy, cz) ;
% 
%     u[0] = atof(cx)/3600*PI/180. ;
%     u[1] = atof(cy)/3600*PI/180. ;
%     u[2] = atof(cz)/3600*PI/180. ;
%     }
% 
%   else { fprintf(stdout, " No gyro measurement (read_gyro).\n") ;
%          return(1) ; }
% 
%   for(i=0;i<3;i++) w[i] = u[i] - b[i] ;
%   return (0) ;
% }
