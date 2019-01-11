function WriteArcGrid(x,y,M,filename)

% WriteArcGrid(x,y,M,'filename') writes the values in matrix M with x and y
% coordinates x and y to the ArcInfo floating point ASCII grid filename.asc. 
% Note that filename must be a character string. See ArcInfo documentation 
% for details on ASCII grid file format. If filename.asc exists in the 
% current directory, it is overwritten.
%
% It is assumed that the elevations are evenly spaced in both directions,
% and that the x and y spacing are the same (the only options in Arc).
%
% It is assumed that x and y have the same ordering as the columns and
% rows, respectively, in M (i.e., the first element in x is the x
% coordinate of the first column in M, and the first element in y is the y
% coordinate of the first row in M.
%
% It is assumed that x is positive east, and y is positive north.

% Copyright (C) 2004-2010 Taylor Perron <perron@mit.edu>
% 
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation. You should have received a copy of the GNU 
% General Public License along with this program.  If not, see 
% http://www.gnu.org/licenses.

if nargin ~= 4, help(mfilename), return, end

% strips the '.asc' from filename if it's there; we'll add it back on later
[sizefnr sizefnc] = size(filename);
if strcmp(filename(sizefnc-3:sizefnc),'.asc')
    filename = filename(1:sizefnc-4);    
end

% if filename.asc exists, delete it
if exist([filename '.asc'],'file')
    delete([filename '.asc'])
end

% write the header file
[ny nx] = size(M);
cellsize = x(2) - x(1);
xllcorner = x(1) - cellsize/2;
yllcorner = y(end) - cellsize/2;


fid=fopen([filename '.asc'],'w');

fprintf(fid, 'ncols %d\nnrows %d\nxllcorner %.6f\nyllcorner %.6f\ncellsize %.6f\nNODATA_value -9999\n', nx, ny, ...
    xllcorner, yllcorner, cellsize);

% Replace NaNs with ArcInfo NODATA value (-9999)
M(isnan(M))=-9999;

% Write the matrix values to the file

% show a progress bar
h = waitbar(0,'Writing ArcInfo grid...');

for i=1:ny
    fprintf(fid, ' %g',M(i,:));
    fprintf(fid, '\n');
    f=i/ny;
    if ~rem(round(f*100),10)
        waitbar(f)
    end
end

close(h);
fclose(fid);
