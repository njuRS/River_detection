function [H Wss] = Hann2D(M)

% [H Wss] = Hann2D(M) 
% 
% Windows matrix M with an elliptical Hann (raised cosine) window. Returns
% the windowed data in H. Also returns the summed square of weighting 
% coefficients, Wss, used in the normalization of the power spectrum.

% Copyright (C) 2004-2010 Taylor Perron <perron@mit.edu>
% 
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation. You should have received a copy of the GNU 
% General Public License along with this program.  If not, see 
% http://www.gnu.org/licenses.

[ny nx] = size(M);
a = (nx+1)/2; b = (ny+1)/2; % matrix coordinates of centroid of M
[X Y] = meshgrid(1:nx,1:ny);

theta = (X==a).*(pi/2) + (X~=a).*atan2((Y-b),(X-a)); % angular polar coordinate

r = sqrt((Y-b).^2 + (X-a).^2); % radial polar coordinate
rprime = sqrt((a^2)*(b^2)*(b^2*(cos(theta)).^2 + a^2*(sin(theta)).^2).^(-1)); % 'radius' of ellipse for this theta

hanncoeff = (r < rprime).*(0.5*(1 + cos(pi*r./rprime)));
H = M.*hanncoeff;

Wss = sum(sum(hanncoeff.^2));
