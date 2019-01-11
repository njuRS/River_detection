function [M F] = SpecFilt2D(M, dx, dy, f, filttype)

% M = SpecFilt2D(M, dx, dy, f, filttype)
%
% Filters a matrix M in the spectral domain.
%
% Input arguments:
%
%     dx, dy   - cell size
%
%     f        - Array of transition frequencies that define the 
%                2D filter.
%
%                For a lowpass filter, f = [flo fhi], where flo is the 
%                frequency at which the filter begins to taper down from 1, 
%                and fhi is the frequency at which the filter tapers close 
%                to zero. 
%                
%                For a bandpass filter, f = [flo1 flo2 fhi1 fhi2], where
%                flo1 is the frequency at which the bandpass filter starts 
%                to increase appreciably above zero, flo2 is the frequency 
%                at which it reaches 1, fhi1 is the frequency at which the 
%                bandpass filter then starts to taper from 1 down to zero, 
%                and fhi2 is the frequency by which it is very nearly zero.
%
%                For a highpass filter, f = [flo fhi], where flo is the 
%                frequency at which the filter starts to increase 
%                appreciably above zero, and fhi is the frequency at which 
%                it reaches 1.
%
%     filttype - String specifying the filter type. Must be 'lowpass', 
%                'highpass', or 'bandpass'.
%
% Output argument:
%
%     M        - The filtered matrix
%     F        - The filter matrix
%
% Dependencies: detrend.m-->lsplane.m, specfilt2d.m

% Copyright (C) 2007-2010 Taylor Perron <perron@mit.edu>
% 
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation. You should have received a copy of the GNU 
% General Public License along with this program.  If not, see 
% http://www.gnu.org/licenses.

[ny nx]=size(M); % number of rows and columns

% get frequency matrix using fftdem
pad=0; window=0;
[Pmat fmat] = fft2D(M, dx, dy, pad, window);

% Remove first-order (planar) trend, but do not window
% M = Detrend(M);

% Do a 2D FFT, padding with zeros.
M = fftshift(fft2(M));

% make the filter matrix
F = Make2DFilt(fmat, f, filttype);

% take the inverse FFT of the filtered spectrum
M = real(ifft2(ifftshift(M.*F)));
M = M(1:ny,1:nx);
