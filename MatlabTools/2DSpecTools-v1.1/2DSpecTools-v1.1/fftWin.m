function [ic jc w q v] = fftWin(M,dy,dx,ywin,xwin,yshift,xshift)

% [ic jc w q v] = fftWin(M,ywin,xwin,yshift,xshift)
%
% Performs a windowed spectral analysis of the matrix M, using a window of
% size [ywin xwin] points. The window is moved in increments of yshift and
% xshift points in the y and x directions.
%
% Returns:
% ic and jc, vectors of row and column indices of the window centers
% w, matrix of weighted peak wavelengths
% q, matrix of orientations of peak wavelengths (degrees CCW from E)
% v, matrix of variances of M within each window

% Copyright (C) 2010 Taylor Perron <perron@mit.edu>
% 
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation. You should have received a copy of the GNU 
% General Public License along with this program.  If not, see 
% http://www.gnu.org/licenses.

% Spectral analysis preferences
pad = 1;
window = 1;

[Ny Nx] = size(M);

% Round xwin and ywin to nearest odd integer
xwin = fix(xwin); ywin = fix(ywin);
if rem(xwin,2)==0, xwin=xwin+1; end
if rem(ywin,2)==0, ywin=ywin+1; end

% Make x and y index vectors that define the window
xwin = (1:xwin)-(xwin+1)/2; 
ywin = (1:ywin)-(ywin+1)/2; 


% Identify window center locations
xshift = round(xshift);
yshift = round(yshift);

ic = 1:yshift:Ny;
jc = 1:xshift:Nx;

% Allocate memory
w = zeros(length(ic),length(jc));
q = zeros(size(w));
v = zeros(size(w));

% Show a progress bar
h = waitbar(0,'Performing windowed spectral analysis...');

for i = 1:length(ic)
    for j = 1:length(jc)
        
        % Discard the fraction of the window that goes off the grid edges
        Z = M(ic(i)+ywin((ic(i)+ywin)>0 & (ic(i)+ywin)<=Ny),...
              jc(j)+xwin((jc(j)+xwin)>0 & (jc(j)+xwin)<=Nx));
        
        % Variance within the window
        v(i,j) = var(Z(:)); 
        
        % 2D FFT
        [P f] = fft2D(Detrend(Z),dx,dy,pad,window);
        
        % Find the highest peak in the spectrum
        [peak, peaki] = max(P(:));

        % Peak wavelength
        %  w(i,j) = 1/f(peaki);

        % Wavelength weighted by power:
        L = 1./f;
        L(f==0)=0; % set DC component to zero rather than Inf
        w(i,j) = sum(sum(L .* P))/sum(P(:));
        
        [Ly Lx] = size(P);
        [ip jp] = ind2sub([Ly Lx],peaki);
        
        % 2D wavenumbers
        [Kx Ky] = meshgrid((-Lx/2):(Lx/2-1),(Ly/2):-1:(-Ly/2+1)); 
        
        % Matrix of orientations
        Q = atan2(Ky*dy,Kx*dx);
        
        % Peak orientation        
        q(i,j) = rad2deg(Q(ip,jp));
        
    end
    
    % Increment progress bar
    waitbar(i/length(ic))
    
end

% Close progress bar
close(h);