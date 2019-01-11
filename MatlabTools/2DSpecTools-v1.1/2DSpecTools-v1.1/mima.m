% mima.m
%
% An exercise prepared for the NSF-sponsored workshop:
% New Tools in Process-Based Analysis of Lidar Topographic Data
% Boulder, Colorado
% June 1, 2010
%
% Figure numbers refer to the accompanying Powerpoint file.
%
% Copyright (c) 2010 Taylor Perron

%% Introduction

% In this exercise, we will analyze a ~1 km^2 portion of a LiDAR dataset
% collected near Merced, CA, by the National Center for Airborne Laser
% Mapping (NCALM). The PIs on the project are Sarah Reed and Ron Amundson
% of UC Berkeley. The main features of interest in the dataset are the
% ubiquitous "mima mounds," (Fig. 1) which are associated with the 
% formation of vernal pools. The terrain contains no trees, and so we will 
% be working with gridded but unfiltered data.
%
% Questions about the mima mounds that we will attempt to answer:
% 
% o What is the characteristic mound spacing?
% o Are there spatial trends in mound size, shape, or organization? 
% o Do these trends correlate with the background topography?
% o How much material has been transported to form the mounds, and how much
%   work was required?
% 
% To address these and other questions, it would be useful to isolate the
% mounds from the background topography. The purpose of this exercise is to
% demonstrate how to use spectral analysis and spectral filtering to
% isolate the mounds and measure their topographic characteristics.
% Specific procedures we will use include:
% 
% o Detrending and windowing of DEMs prior to performing spectral analysis
% o Use of the 2D Fast Fourier Transform to estimate 2D power spectra
% o Identification of quasiperiodic topographic features and measurement of
%   their characteristics, including wavelength, amplitude, and orientation 
% o Spectral-domain filtering to isolate features of interest
% o Windowed spectral analysis to map spatial variations in landform
%   characteristics
% o Optimal filtering to remove noise


%% %%%%%%%%%%%%
% 0. SETTINGS %
%%%%%%%%%%%%%%%

%% A note on cell mode

% This script uses Matlab's cell mode, which is why the section where your
% cursor currently resides is highlighted. The easiest way to progress
% through this tutorial is with the following commands, which are also
% available under the "Cell" menu:
%
% Cmd + Enter: Evaluate all the commands in the current cell
% or
% Shift + Cmd + Enter: Evaluate current cell and advance to the next cell
%
% (On a PC, replace Cmd with Ctrl)

%% Set up the Matlab desktop

% First, take a moment to dock the Editor so you can see this script and
% the command line at the same time.

% Make new figures docked and white
set(0,'DefaultFigureWindowStyle','docked','DefaultFigureColor','w')


%% %%%%%%%%%%%%%
% 1. LOAD DATA %
%%%%%%%%%%%%%%%%

%% Import gridded data

% The LiDAR data are stored as an ArcGIS-formatted ASCII grid (*.asc)
[Z, dim] = ReadArcGrid('sub_orthowv01_20150718160144_for_wq7_new.asc');

% The struct variable dim contains fields describing the coordinates and
% dimensions of Z. These are based on the information extracted from the
% grid file header.
dx = abs(dim.x(2) - dim.x(1)); % grid spacing in the x-direction
dy = abs(dim.y(2) - dim.y(1)); % grid spacing in the y-direction
[Ny Nx] = size(Z); % grid dimensions

%% View data

% Display a shaded relief map (Fig. 2)
figure('Name','Fig. 2: Shaded relief','NumberTitle','off')
ShadePlot(dim.x,dim.y,Z)

% Use the magnifying glass tool to zoom in for a closer look at the mounds.
% What is the typical inter-mound spacing? Note that the grid spacing is 
% 1m. They're about 10 pixels apart, so we'll be on the lookout for a
% feature at a wavelength of ~10m in the power spectrum.


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. CALCULATE POWER SPECTRUM %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preprocessing

% Before calculating the power spectrum, we'll make some decisions about
% the preprocessing steps we'll take. If the DEM has a non-zero mean or a
% background slope across the entire grid, this will contaminate the
% spectrum with long-wavelength signals. So our first step will be to 
% detrend the DEM by fitting a least-squares plane to the elevations and 
% then subtracting this fit.
Zo = Z; % Save the original elevations for later
Z = Detrend(Z);
plane = Zo - Z; % Save the least-squares plane for re-trending later

% Second, the fast Fourier transform proceeds fastest if the dimensions of 
% the input matrix are integer powers of two, which we can achieve by 
% padding the DEM with zeros. 
pad = 1; % 1 means pad the data, 0 no padding.

% Third, because the edges of our DEM are not perfectly periodic, the
% spectrum can become contaminated by frequencies used to "fit" the edge
% discontinuity. We can mitigate this effect by multiplying the DEM by a
% function that tapers to zero at the edges.
window = 1; % 1 means window the data, 0 no window

%% 2D FFT

% Calculate the power spectrum using the 2D Fast Fourier Transform (FFT).
% Note that the zero padding and windowing happen inside fft2D.
[Pm, fm Pv fv] = fft2D(Z,dx,dy,pad,window); 

% A few words about the output from the previous step:
% 
% We have calculated the Discrete Fourier Transform (DFT) periodogram. This
% is different from the power spectral density, which is also commonly used
% as an estimate of the power spectrum. 
%
% Pm is the 2D power spectrum matrix, and fm is the corresponding matrix of 
% radial frequencies. Pv is a vector version of the power spectrum, and fv 
% is the corresponding vector of frequencies. Units of Pm and Pv are 
% [units of Z]^2, which in our case is length^2 since Z is a matrix of 
% elevations. Units of fm and fv are [units of x]^-1, which for us is 
% 1/length since x is distance. Figs. 3a-e show examples of a simple 1D
% signal and its power spectrum, and Fig. 3f shows an example of a simple
% surface and its 2D power spectrum.
%
% What do the units mean? The periodogram is a measure of how much of the
% original elevation field's variance falls within a given frequency range.
% You can check that the sum of the periodogram is roughly equal to the
% variance in Z. (It will be somewhat less due to the zero padding.)
%
% What about the frequency limits? Wavelength is equal to the inverse of 
% frequency. The shortest wavelength that can be resolved by 2D data is 
% oriented diagonally, and is equal to sqrt(dx^2 + dy^2), or sqrt(2) for 
% our data (Fig. 4). The longest wavelength that can be resolved is 
% infinite in principle, but in practice we wouldn't feel comfortable 
% trying to detect any signal with a wavelength longer than our dataset. To 
% be even more conservative, the longest wavelength we trust should be a 
% few times shorter than the width of the DEM -- a few hundred meters in 
% this case.

% Clean up workspace
clear pad window

%% %%%%%%%%%%%%%%%%%%%%%%%%
% 3. ANALYZE THE SPECTRUM %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1D Power Spectrum

% Plot the 1D version of the spectrum (Fig. 5). We'll plot the 2D spectrum 
% a few steps later, but for now it's easier to visualize in 1D.
s1d = figure('Name','Fig. 5: 1D Spectrum','NumberTitle','off');
subplot(2,1,1)
[axf axw] = SpecPlot1D(fv,Pv);

% Do we see anything at a 10m wavelength? Yes, there is a broad peak. The 
% shorter-wavelength peaks are just pixel-scale noise.
hold(axw,'on')
plot([10 10],get(axw,'ylim'),'k')


%% Background Spectrum

% Now let's see what the spectrum looks like in 2 dimensions. If we just
% looked at it as is, we wouldn't see much: as the 1D spectrum shows,
% longer-wavelength signals have much higher power, and would swamp any
% shorter-wavelength detail. To make the 2D spectrum easier to view, we'll
% remove the power-law background trend. Because there are so many more
% points at higher frequencies, we'll bin the 1D spectrum and fit the trend
% to the binned values (Fig. 5).

nbin = 20; % Number of bins
B = bin(log10(fv),log10(Pv),nbin,0); % Bin the log-transformed data

% Plot the binned values
hold(axf,'on')
plot(axf,10.^B(:,1),10.^B(:,2),'ok','markerfacecolor','w')

% Fit a trend with the form P ~ 1/f^n, and plot it
fit = robustfit(B(:,1),B(:,2));
plot(axf,10.^B(:,1),10^fit(1)*(10.^B(:,1)).^fit(2),'k')

%% 2D Power Spectrum

% Use this fit to normalize the 2D spectrum
Pmn = Pm./(10^fit(1)*fm.^fit(2));

% Plot the normalized 2D spectrum (Fig. 6)
figure('Name','Fig. 6: 2D Spectrum','NumberTitle','off')
SpecPlot2D(fm,log10(Pmn));

% The donut of elevated values at radial frequencies of ~0.1 is the
% spectral signature of the mima mounds. It's a donut because the mounds
% don't have a strongly preferred orientation over the entire landscape. 
% The two peaks very close to the center of the square correspond to the 
% NW-SE-trending valleys, and the sharp features aligned with these peaks 
% at higher frequencies are probably harmonics.

% Clean up workspace
clear Pm fm Pv nbin B axf axw fit Pmn

%% %%%%%%%%%%%%%
% 4. FILTERING %
%%%%%%%%%%%%%%%%

%% High-pass filter

% Now let's isolate the signal associated with the mima mounds. From our 
% 1D and 2D spectra, it appears that the peak associated with the mounds
% begins at frequencies of ~0.03 to 0.05 m^-1 (wavelengths of 20 to 30m). 
% The mounds are only a few pixels wide, so we want to retain all 
% frequencies higher than about 0.04. We can accomplish this with a
% high-pass filter.

% Define the transition frequencies for the filter: floHP is the frequency 
% at which the high-pass filter starts to increase appreciably above zero. 
% fhiHP is the frequency at which it reaches 1.
%floHP = 1/100; fhiHP = 1/15;
floHP = 1/100; fhiHP = 1/2;

% Plot the filter alongside the 1D spectrum to illustrate its shape (Fig. 
% 5). Note how the filter transitions smoothly from 0 to 1 -- the shape is 
% half a Gaussian in linear frequency.
figure(s1d)
subplot(2,1,2)
semilogx(fv,Make2DFilt(fv, [floHP fhiHP], 'highpass'),'b')
set(gca,'box','off','tickdir','out','ylim',[0 1.01])
xlabel('Radial frequency (m^{-1})')
ylabel('Filter')

% Filter the DEM
Zhp = SpecFilt2D(Z,dx,dy,[floHP fhiHP],'highpass');

%% View filtered DEM

% Display a shaded relief map of the high-pass-filtered data (Fig. 7)
hp = figure('Name','Fig. 7: High-pass-filtered DEM','NumberTitle','off');
ShadePlot(dim.x,dim.y,Zhp)

%% Low-pass filter

% What's left over is the low-pass-filtered topography (Fig. 8)
Zlp = Z - Zhp;
lp = figure('Name','Fig. 8: Low-pass-filtered DEM','NumberTitle','off');
ShadePlot(dim.x,dim.y,Zlp)
% We could also do the low-pass filtering with the following:
% floLP = floHP; fhiLP = fhiHP;
% Zlp = SpecFilt2D(Z,dx,dy,[floLP fhiLP],'lowpass');


%% Remove channels

% Note that the mima mounds are not the only features in the
% high-pass-filtered topography: channels are also short-wavelength
% features. We could try to use a band-pass filter to eliminate the
% channels, but their width is close enough to the mound wavelength that we
% would risk eliminating important details of the mounds. Instead, we'll
% try to remove them in the spatial domain by exploiting the fact that they
% are concave-up features with low elevations. 

% Before calculating derivatives, "re-trend" the data by adding the 
% least-squares plane back in
Zlp = Zlp + plane;

% Calculate derivatives of the low-pass-filtered topography 
[DX DY] = gradient(Zlp,dx,dy); 
DX = -DX; % Correct for the fact that y is positive up
G = sqrt(DX.^2 + DY.^2); % Magnitude of the gradient
D = atan2(DY,DX); % Direction of the gradient vector (aka aspect)
C = 4*del2(Zlp,dx,dy); % Laplacian

% Create a binary grid of channel locations and set elevations in those
% locations equal to the mean of the rest of the grid.
Cmax = 0.01;
Zmax = 0;
chan = C > Cmax & Z < Zmax;

%Zhp(chan) = mean(Zhp(~chan)); %beautiful!  DON'T DELETE CHANNELS 10.2.2015

% Display a shaded relief map of the result (Fig. 7)
figure(hp)
ShadePlot(dim.x,dim.y,Zhp)

%% Remove edge effects

% You can also see some edge effects produced by the removal of long 
% wavelengths. Let's clip these off.
c = 30;
xclip = c:Nx-c; yclip = c:Ny-c;
Zhp = Zhp(yclip,xclip);
Zlp = Zlp(yclip,xclip);
xhp = dim.x(xclip);
yhp = dim.y(yclip);

% Display the result (Figs. 7 & 8)
figure(hp)
ShadePlot(xhp,yhp,Zhp)
figure(lp)
ShadePlot(xhp,yhp,Zlp)

% Do the same to the other grids we'll use below
chan = chan(yclip,xclip);
G = G(yclip,xclip);
D = D(yclip,xclip);

% Clean up workspace
clear fv floHP fhiHP s1d Cmax Zmax DX DY c xclip yclip hp lp


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. WINDOWED SPECTRAL ANALYSIS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Windowed 2D FFT

% Now that we have isolated the mounds from the background topography, we
% can investigate how they vary across the landscape. The spectral analysis
% we performed previously only provides a measure of average mound
% characteristics across the entire DEM, but the lidar data is sufficiently
% dense that we can compute multiple spectra within smaller windows across
% the DEM, and derive a measure of mound characteristics within each window.
xwin = 100; ywin = 100; % Window dimensions in pixels
xstep = 10; ystep = 10; % Distances between successive window positions
[ic jc wavelength angle variance] = ...
    fftWin(Zhp,dy,dx,xwin,ywin,xstep,ystep);


% Visualize the results of the windowed analysis by overlaying them on
% a shaded relief map of the original DEM:
opacity = 0.25;

%% Variance

% Variance provides a local measure of mound amplitude (Fig. 9):
figure('Name','Fig. 9: Variance map','NumberTitle','off');
CShadePlot(dim.x,dim.y,Z,xhp(jc),yhp(ic),variance,opacity); 
title('Variance (m^2)')

%% Wavelength

% Spatial variability in the mound wavelength (Fig. 10):
figure('Name','Fig. 10: Wavelength map','NumberTitle','off');
wmax = 15;
wavelength(wavelength>wmax) = wmax; % Avoid scale saturation in channels
CShadePlot(dim.x,dim.y,Z,xhp(jc),yhp(ic),wavelength,opacity); 
title('Mound wavelength (m)')

%% Orientation

% Direction of mound elongation (Fig. 11):
angle = angle+90; % Convert wave direction to mound elongation direction
angle(angle<0) = angle(angle<0)+360;
angle = angle - 180; % Make all values 0 to 180
figure('Name','Fig. 11: Orientation map','NumberTitle','off');
CShadePlot(dim.x,dim.y,Z,xhp(jc),yhp(ic),angle,opacity); 
title('Mound orientation (degrees CCW from E)')

% Clean up workspace
clear xwin ywin xstep ystep wmax


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6. ANALYSIS OF MOUND CHARACTERISTICS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Does mound height depend on slope?
G = G(ic,jc);
figure('Name','Fig. 12: Mound height vs. Slope','NumberTitle','off');
chan = chan(ic,jc); % A channel map at the windowed resolution

% There are so many lidar points that it's useful to include a point
% density map to make trends clearer:
scattercloud(G(~chan),variance(~chan),25,1,'.k','jet');
xlabel('Gradient')
ylabel('Variance (m^2)');
set(gca,'tickdir','out')
% We can see that mound height is typically highest on gentle slopes, and
% declines as slopes steepen (Fig. 12). 

%% Does mound height depend on aspect?
D = D(ic,jc);
figure('Name','Fig. 13: Mound height vs. Aspect','NumberTitle','off');
h1 = polar(D(~chan),variance(~chan),'.');
set(h1,'MarkerEdgeColor',[.5 .5 .5])
d = bin(D(:),variance(:),24,0);
hold on
h2 = polar(d(:,1),d(:,8),'-ok');
set(h2,'MarkerFaceColor','y','MarkerSize',10)
title('Mound variance (m^2) vs. slope aspect (degrees CCW from E)')
% Mounds are taller on NE-facing slopes, and shorter on SW-facing slopes
% (Fig. 13).

%% We can make the trend more apparent by looking only at the binned values:
set(h1,'visible','off')


%% Are mounds aligned or elongated in the downslope direction? 
O = D;  
O(O<0)=O(O<0)+pi; % Make all slope orientations 0 - 180 deg
Gmin = 0.1; % Minimum gradient -- we don't expect a strong relationship on
            % gentle slopes
figure('Name','Fig. 14: Mound orientation vs. Aspect','NumberTitle','off');
scattercloud(rad2deg(O(~chan & G>Gmin)),angle(~chan & G>Gmin),25,1,'.k','jet');
axis square
hold on
plot([0 180],[0 180],'w')
xlabel('Slope orientation (deg CCW from E)');
ylabel('Mound orientation (deg CCW from E)');
title('Gradient > 0.1')
set(gca,'tickdir','out','xlim',[0 180],'ylim',[0 180])
% It's noisy, but the clustering around the 1:1 line suggests that the
% mounds are preferentially elongated in the downslope direction (Fig. 14).



%% How much energy is required to form a mound? 
% If mounds are formed by transporting soil from troughs to peaks, then a 
% rough estimate of the work W required to build a mound is
% W = rho*g*V*h
% where rho is soil density, g is gravity, V is the volume of soil moved
% from trough to peak, and h is the vertical distance the soil's center of
% mass was displaced, which is roughly equal to half the mound height:
%         __
%  /\  /\ __ h
% /  \/  \
%
% If the mound wavelength is lambda, then V ~ h*lambda^2. So our estimate
% of the work required to build a mound is:
rho = 2200; % kg m^-3
g = 10; % m s^-2
h = sqrt(variance); % Standard deviation (m), an estimate of h
W = rho*g*h.^2.*wavelength.^2; % J

% Plot energy expenditure over shaded relief (Fig. 15)
figure('Name','Fig. 15: Energy expenditure','NumberTitle','off');
CShadePlot(dim.x,dim.y,Z,xhp(jc),yhp(ic),W/1000,opacity); 
title('Energy expenditure (kJ)')

% Clean up workspace
clear G C D O chan s d h1 h2 Gmin rho g opacity


%% %%%%%%%%%%%%%%
% 7. DE-NOISING %
%%%%%%%%%%%%%%%%%

%% A noisy DEM

% What if we want to eliminate, rather than isolate, a component of a
% surface? Noise removal is a common filtering operation, and there is an
% optimal way to construct the filter. 

% Add some Gaussian noise to the elevations in the spatial domain to create 
% a noisy DEM (Fig. 16)
sigma = 0.2; % Standard deviation of the noise in meters
No = Zo + sigma*randn(Ny,Nx); % Noisy elevations
figure('Name','Fig. 16: Noisy DEM','NumberTitle','off'); 
ShadePlot(dim.x,dim.y,No);

%% Power spectrum

% Get the 2D Power spectrum using the same steps outlined above
N = Detrend(No);
plane = No - N; % Save the least-squares plane for re-trending later
pad = 1;
window = 1;
[Pm, fm Pv fv] = fft2D(N,dx,dy,pad,window); 

% Plot the 1D version of the spectrum (Fig. 17)
figure('Name','Fig. 17: 1D Spectrum of noisy DEM','NumberTitle','off');
axf = SpecPlot1D(fv,Pv);

nbin = 20;
B = bin(log10(fv),log10(Pv),nbin,0);

hold(axf,'on')
plot(axf,10.^B(:,1),10.^B(:,2),'ok','markerfacecolor','w')


%% Noise estimation

% Next, we manually identify parts of the spectrum that appear to be
% dominated by noise. The high-frequency end of the spectrum is flat, the
% expected trend for white noise (Fig. 17).
imeas = 1:20; % Bins that correspond to the measured spectrum (all of them)
inoise = 17:20; % Bins that are dominated by noise

% Extrapolate the inferred noise spectrum to parts of the spectrum that are
% dominated by signal. In this case, we can capture the assumed shape of 
% the noise spectrum pretty well with a constant value.
noise = 10.^min(B(inoise,2)) * ones(size(fm));

% Fit a smooth surface to the measured spectrum by interpolating linearly
% between the binned values.
measured = interp1(B(imeas,1),B(imeas,2),log10(fm));
measured = 10.^measured; % Transform back to linear spectral power
measured(isnan(measured)) = 0; % Elements that couldn't be interpolated


% Plot the measured and inferred spectra (Fig. 17)
plot(axf,fm(:),measured(:),'+k')
plot(axf,fm(:),noise(:),'+b')

%% Filter construction

% Estimate the true signal spectrum by subtracting the estimated noise
% spectrum from the measured spectrum. The optimal filter is then 
% constructed by taking the ratio of the estimated signal spectrum to the
% measured spectrum. The filter will be nearly 1 at frequencies where the 
% spectrum is dominated by signal, and nearly zero where it is dominated by
% noise -- just what we want.
signal = measured - noise; 
F = signal./measured;
F(F==Inf | F==-Inf) = 0; % Correct elements where we just divided by zero

%% Noise removal

% Now we do the filtering. This time we don't window the DEM, because we
% want to preserve the actual elevations. 
padsize = 1024; % Next highest power of 2
N = fftshift(fft2(N,padsize,padsize)); % Fourier transform of noisy DEM
N = ifft2(ifftshift(N.*F)); % Inverse transform of filtered FT
N = plane + N(1:Ny,1:Nx); % Clip off the padding from the de-noised DEM and
                          % re-trend by adding the least-squares plane

%% View de-noised DEM

% Plot the de-noised DEM (Fig. 18)
figure('Name','Fig. 18: De-noised DEM','NumberTitle','off');
ShadePlot(dim.x,dim.y,N);

% Clean up workspace
clear sigma plane pad window Pm fm Pv fv axf nbin B imeas inoise noise measured signal F padsize 


%% %%%%%%%%%%%%%%%
% 8. SAVE OUTPUT %
%%%%%%%%%%%%%%%%%%

%% Save workspace
save mima

%% Export grids

% Should we want to use any of the filtered grids in a GIS, we can use the
% included function that writes an Arc-formatted ASCII grid. This grid
% can then be imported directly into ArcMap, or converted to ESRI's binary
% grid format using the ASCIIGRID command in ArcInfo or the "ASCII to
% Raster" function in ArcToolbox. Examples:

% The high-pass-filtered topography
WriteArcGrid(xhp,yhp,Zhp,'merced_highpass')

% The map of mound wavelength
WriteArcGrid(xhp(jc),yhp(ic),wavelength,'wavelength')