function ShadePlot(x,y,Z)

% ShadePlot(x,y,Z)
%
% Plots a shaded relief map of the elevation matrix Z with x and y 
% coordinates x and y.
%
% Dependencies: hillshade.m

% Copyright (C) 2010 Taylor Perron <perron@mit.edu>
% 
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation. You should have received a copy of the GNU 
% General Public License along with this program.  If not, see 
% http://www.gnu.org/licenses.

imagesc(x,y,hillshade(Z,x,y)); 
axis image
colormap gray
set(gca,'ydir','normal','tickdir','out')
