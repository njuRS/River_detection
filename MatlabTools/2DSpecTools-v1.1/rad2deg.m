function angleInDegrees = rad2deg(angleInRadians)
% RAD2DEG Convert angles from radians to degrees
%
%   angleInDegrees = RAD2DEG(angleInRadians) converts angle units from
%   radians to degrees.
%
%   See also: deg2rad, fromDegrees, fromRadians, toDegrees, toRadians.

% Copyright 2007 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2007/02/11 05:47:41 $

angleInDegrees = (180/pi) * angleInRadians;
