function CShadePlot(x,y,Z,xc,yc,C,opacity)

% CShadePlot(x,y,Z,C)
%
% Plots shaded relief of the elevation matrix Z with x and y coordinates x
% and y, and overlays transparent colored plot of matrix C with x and y
% coordinates xc and yc. opacity (0 to 1) sets the opacity of the overlay.
%
% Dependencies: hillshade.m
%
% Method based on Mathworks technical note:
% http://www.mathworks.com/support/tech-notes/1200/1215.shtml

% Copyright (C) 2004-2010 Taylor Perron <perron@mit.edu>
% 
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation. You should have received a copy of the GNU 
% General Public License along with this program.  If not, see 
% http://www.gnu.org/licenses.

H = hillshade(Z,x,y);

colormap([gray(64);jet(64)])

h(1) = image(x,y,H); axis image
hold on
h(2) = image(xc,yc,C); axis image

m = 64;  % 64 elements in each colormap

% CData for hillshade
cminH = min(H(:));
cmaxH = max(H(:));
CH = min(m,round((m-1)*(H-cminH)/(cmaxH-cminH))+1);

% CData for overlay
cminC = min(C(:));
cmaxC = max(C(:));
CC = min(m,round((m-1)*(C-cminC)/(cmaxC-cminC))+1);
CC = 64+CC;

set(h(1),'CData',CH);
set(h(2),'CData',CC);
caxis([min(CH(:)) max(CC(:))])

alpha(h(2),opacity)

set(gca,'ydir','normal','tickdir','out')

hold off


% Create colorbar for overlay
c = C(:);
ntick = 5;
tickloc = (1:(64-1)/ntick:64)+64;
ticklab = (min(c):(max(c)-min(c))/ntick:max(c));
colorbar('tickdir','out','box','off','ylim',[65 128],...
         'ytick',tickloc,'yticklabel',ticklab);
