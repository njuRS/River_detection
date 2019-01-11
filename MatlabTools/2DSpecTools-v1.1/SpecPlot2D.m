function ax = SpecPlot2D(f,P)

% ax = SpecPlot2D(f,P)
%
% Plots a two-dimensional power spectrum with frequency matrix f and power
% spectrum P. Returns a handle to the axes.

% Copyright (C) 2004-2010 Taylor Perron <perron@mit.edu>
% 
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation. You should have received a copy of the GNU 
% General Public License along with this program.  If not, see 
% http://www.gnu.org/licenses.

ax = axes;
imagesc(P); axis image
xlabel('x frequency')
ylabel('y frequency')

[nfy nfx] = size(f);
nyq = f(nfy/2+1,1); % the Nyquist frequency in the x direction

% Note that these frequency labels are approximate due to the offset of the
% DC (zero frequency) element. The frequencies in f are exact.
numticks=5;
set(gca,'XTick',linspace(1,nfx,numticks))
set(gca,'YTick',linspace(1,nfy,numticks))
set(gca,'XTicklabel',linspace(-nyq,nyq,numticks))
set(gca,'YTicklabel',linspace(nyq,-nyq,numticks))
set(gca,'TickDir','out')
