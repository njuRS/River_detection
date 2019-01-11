function s = bin(x,y,nbin,overlap)

% s = bin(x,y,numwin)
%
% Bins the function y(x) using nbin bins of equal x-width (as opposed to 
% bins containing a constant number of data points, as is the case with 
% many other routines). 
%
% Input arguments:
%    x        - a vector containing the values of the independent variable
%    y        - a vector containing the values of the dependent variable y
%               at the values in x. length(y) == length(x).
%    nbin     - number of bins of equal x-width
%    overlap  - a flag that determines whether the bins will be 
%               non-overlapping (overlap=0, default) or 50% overlap 
%               (overlap=1)
%
% Output arguments:
%    s        - matrix consisting of nbin rows and 8 columns: 
%               [x (center of bin), mean, standard deviation, standard 
%               error, n, max, min, median]

% Copyright (C) 2004-2010 Taylor Perron <perron@mit.edu>
% 
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation. You should have received a copy of the GNU 
% General Public License along with this program.  If not, see 
% http://www.gnu.org/licenses.

if nargin<4, overlap=0; end

% if x and y are row vectors, make them column vectors
x = x(:);
y = y(:);

% sort x and y by x
sorted = sortrows([x y],1);
x = sorted(:,1); y = sorted(:,2);

% find the extrema of x
xmin = x(1); xmax = x(end);
xrange = xmax - xmin;


% determine the window width
if overlap,
    w = 2*xrange/(nbin+1); % for windows with 50% overlap (smoother)
else
    w = xrange/nbin; % for windows with no overlap
end

% Allocate memory for the binned data
s = zeros(nbin,8);

% loop through the bins
for i=1:nbin
	
    % determine min and max x values of current window position
    if overlap
        xlo = xmin+(i-1)*(w/2); % for windows with 50% overlap
    else
        xlo = xmin+(i-1)*w; % for windows with no overlap
    end
    xhi = xlo+w;
	
    % find min and max indices of x vector corresponding to this range
	window = find((x >= xlo) & (x <= xhi));
    mini = min(window); maxi = max(window);
    
    % calculate mean, standard dev, standard error, and n of points that 
    % fall within this window, but watch out for windows with only one 
    % point:
    if isempty(window)
        s(i,:) = [mean([xlo xhi]) 0 0 0 0 0 0 0];
    else
        s(i,:) = [mean([xlo xhi]) mean(y(mini:maxi)) ...
                 std(y(mini:maxi)) std(y(mini:maxi))/sqrt(maxi-mini+1)...
                 maxi-mini+1 max(y(mini:maxi)) min(y(mini:maxi))...
                 median(y(mini:maxi))];
    end
    
end
