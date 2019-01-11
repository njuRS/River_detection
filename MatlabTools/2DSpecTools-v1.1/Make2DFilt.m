function F = Make2DFilt(fmat, f, filttype)

% F = Make2DFilt(fmat, f, filttype)
%
% Constructs a 2D spectral filter. See SpecFilt2D.m for explanation.

% Copyright (C) 2007-2010 Taylor Perron <perron@mit.edu>
% 
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation. You should have received a copy of the GNU 
% General Public License along with this program.  If not, see 
% http://www.gnu.org/licenses.


if nargin < 3
    help(mfilename), 
    return
end

switch filttype
    
    case 'lowpass'
        flo = f(1); fhi = f(2);
        mu=flo;
        sigma=abs(fhi-flo)/3;
        F=Gaussian(fmat,mu,sigma);
        F(fmat<flo)=1;        

    case 'highpass'
        flo = f(1); fhi = f(2);
        mu=fhi;
        sigma=abs(fhi-flo)/3;
        F=Gaussian(fmat,mu,sigma);
        F(fmat>=fhi)=1;

    case 'bandpass'
        flo1 = f(1); flo2 = f(2);
        fhi1 = f(3); fhi2 = f(4);        
        sigmalo = abs(flo2-flo1)/3;
        sigmahi = abs(fhi2-fhi1)/3;
        mulo=flo2;
        muhi=fhi1;
        Flo=Gaussian(fmat,mulo,sigmalo);
        Fhi=Gaussian(fmat,muhi,sigmahi);
        F = Flo.*(fmat<=mulo) + Fhi.*(fmat>=muhi) + 1*(fmat>mulo & fmat<muhi);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function G = Gaussian(freqmat,mu,sigma)

G=exp(-(freqmat-mu).^2/(2*sigma^2));
G=G/max(G(:));