function gb=gaborfilter(sigma,theta,lambda,psi,elongation,type)
 
sigma_x = sigma;
sigma_y = sigma * elongation;
 
% Bounding box
nstds = 3;
xmax = max(abs(nstds*sigma_x*cos(theta)),abs(nstds*sigma_y*sin(theta)));
xmax = ceil(max(1,xmax));
ymax = max(abs(nstds*sigma_x*sin(theta)),abs(nstds*sigma_y*cos(theta)));
ymax = ceil(max(1,ymax));
xmin = -xmax; ymin = -ymax;
[x,y] = meshgrid(xmin:xmax,ymin:ymax);
 
% Rotation 
x_theta=x*cos(theta)+y*sin(theta);
y_theta=-x*sin(theta)+y*cos(theta);

if strcmp(type, 'even')
    gb= 1/(2*pi*sigma_x *sigma_y) * exp(-.5*(x_theta.^2/sigma_x^2+y_theta.^2/sigma_y^2)).*cos(2*pi/lambda*x_theta+psi);
else
    gb= 1/(2*pi*sigma_x *sigma_y) * exp(-.5*(x_theta.^2/sigma_x^2+y_theta.^2/sigma_y^2)).*sin(2*pi/lambda*x_theta+psi);
end


