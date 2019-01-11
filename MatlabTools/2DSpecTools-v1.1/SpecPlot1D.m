function [axf axw] = SpecPlot1D(f,P)

% [axf axw] = SpecPlot1D(f,P)
%
% Plots spectral power or periodogram P against frequency f and wavelength
% 1/f using two x-axes. Returns handles to the frequency and wavelength
% axes.

% Copyright (C) 2004-2010 Taylor Perron <perron@mit.edu>
% 
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation. You should have received a copy of the GNU 
% General Public License along with this program.  If not, see 
% http://www.gnu.org/licenses.

loglog(f,P,'or','markersize',3);
axf = gca;

frange = get(axf,'xlim');
wrange = 1./fliplr(frange);

axw=axes('Position',get(axf,'Position'),...
         'XAxisLocation','top',...
         'YAxisLocation','right',...
         'Color','none',...
         'XColor','k','YColor','k',...
         'xlim',wrange,'xdir','reverse',...
         'ytick',[],'xscale','log');

set(axf,'box','off','tickdir','out')
set(axw,'box','off','tickdir','out')

set(get(axf,'xlabel'),'string','Radial frequency (m^{-1})')
set(get(axf,'ylabel'),'string','Mean squared amplitude (m^2)')
set(get(axw,'xlabel'),'string','Wavelength (m)')
