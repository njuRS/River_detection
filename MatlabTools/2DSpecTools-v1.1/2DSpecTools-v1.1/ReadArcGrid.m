function [M, dimensions] = ReadArcGrid(filename)

% [M dimensions] = ReadArcGrid('filename') 
%
% Reads the ArcInfo ASCII grid filename.asc into a matrix of grid values, 
% M, and a struct array, dimensions, with elements: 
%
% dimensions.ncols     (# of columns in grid)
% dimensions.nrows     (# of rows in grid)
% dimensions.x         (x coordinates of centers of pixels)
% dimensions.y         (y coordinates of centers of pixels)
% dimensions.cellsize  (cell size)
% 
%
% Note: assumes that the .asc file header is in the standard Arc format:
%
% NCOLS xxxxxx
% NROWS xxxxxx
% XLLCORNER xxxxxx
% YLLCORNER xxxxxx
% CELLSIZE xxxxxx
% NODATA_VALUE xxxxxx
%
% If the name of the grid is filename.asc, the input argument can be either
% 'filename.asc' or 'filename'. 

% Copyright (C) 2004-2010 Taylor Perron <perron@mit.edu>
% 
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation. You should have received a copy of the GNU 
% General Public License along with this program.  If not, see 
% http://www.gnu.org/licenses.

if (nargin ~= 1), help(mfilename), return, end

[sizefnr sizefnc] = size(filename);
if strcmp(filename(sizefnc-3:sizefnc),'.asc')
    filename = filename(1:sizefnc-4);    
end

if ~exist([filename '.asc'],'file')
    error(['File ' filename '.asc does not exist.']);
end

% open the file and check the first string. Should read "ncols". If it doesn't,
% tell the user it isn't a valid file & bail out.

fid = fopen([filename '.asc'],'r');
if ~strcmpi(fscanf(fid, ' %s',1),'ncols')
    fclose(fid);
    error(['File ' filename '.asc is not an ArcInfo ASCII grid file.']);
end

% read in the grid description and calculate cellsize etc.
% we are now right after 'ncols' in the file
dimensions.ncols = fscanf(fid, '%g', 1);
dimensions.nrows = fscanf(fid, '%*s%g', 1);
xllcorner = fscanf(fid, '%*s%g', 1);
yllcorner = fscanf(fid, '%*s%g', 1);
dimensions.cellsize = fscanf(fid, '%*s%g', 1);
nodata_value = fscanf(fid, '%*s%g', 1);
dimensions.x = xllcorner + dimensions.cellsize * (0.5 + (0:dimensions.ncols-1));
dimensions.y = flipud(yllcorner + dimensions.cellsize * (0.5 + (0:dimensions.nrows-1)'));

% read the grid values into the variable M
M = fscanf(fid,'%g',[dimensions.ncols dimensions.nrows]);
fclose(fid);
M = M';

% Replace nodata values with NaNs
M(M==nodata_value)=NaN;
