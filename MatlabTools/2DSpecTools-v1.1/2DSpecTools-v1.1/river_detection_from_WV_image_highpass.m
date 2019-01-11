% 0. SETTINGS %
set(0,'DefaultFigureWindowStyle','docked','DefaultFigureColor','w')
% 1. LOAD DATA %
[Z, dim] = ReadArcGrid('sub_orthowv02_12jul231608563_band8.asc');
dx = abs(dim.x(2) - dim.x(1)); % grid spacing in the x-direction
dy = abs(dim.y(2) - dim.y(1)); % grid spacing in the y-direction
[Ny Nx] = size(Z); % grid dimensions

if(mod(Ny,2)==1)   %fft2D cannot process odd number, so preprocess is required here; debug: 5/19/2016
    Ny=Ny-1;
end
if(mod(Nx,2)==1)
    Nx=Nx-1;
end
Z=Z(1:Ny,1:Nx);

Z(isnan(Z)==1)=0; %sometimes there are NAN values near image boundaries

%% View data
%figure('Name','Fig. 2: Shaded relief','NumberTitle','off')
%ShadePlot(dim.x,dim.y,Z)


% 2. CALCULATE POWER SPECTRUM %
%% Preprocessing
Zo = Z; % Save the original elevations for later
Z = Detrend(Z);
plane = Zo - Z; % Save the least-squares plane for re-trending later
pad = 1; % 1 means pad the data, 0 no padding.
%% 2D FFT
[Pm, fm Pv fv] = fft2D(Z,dx,dy,pad,window); 
clear pad window



% 3. ANALYZE THE SPECTRUM %
%% 1D Power Spectrum
s1d = figure('Name','Fig. 5: 1D Spectrum','NumberTitle','off');
subplot(2,1,1)
[axf axw] = SpecPlot1D(fv,Pv);
hold(axw,'on')
plot([10 10],get(axw,'ylim'),'k')


%% Background Spectrum
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
% Clean up workspace
clear Pm fm Pv nbin B axf axw fit Pmn


%%
% 4. FILTERING %
% %% High-pass filter
% floHP = 1/100; fhiHP = 1/1;
% %floHP = 1/13; fhiHP = 1/11;
% 
% % Plot the filter alongside the 1D spectrum to illustrate its shape (Fig. 
% % 5). Note how the filter transitions smoothly from 0 to 1 -- the shape is 
% % half a Gaussian in linear frequency.
% figure(s1d)
% subplot(2,1,2)
% semilogx(fv,Make2DFilt(fv, [floHP fhiHP], 'highpass'),'b')
% set(gca,'box','off','tickdir','out','ylim',[0 1.01])
% xlabel('Radial frequency (m^{-1})')
% ylabel('Filter')
% 
% % Filter the DEM
% Zhp = SpecFilt2D(Z,dx,dy,[floHP fhiHP],'highpass');
% % Low-pass filter
% Zlp = Z - Zhp;
% %Zlp = SpecFilt2D(Z,dx,dy,[floLP fhiLP],'lowpass');

%%
%f=[1/200 1/50 1/10 1/1];
f=[1/1000 1/10];
figure(s1d)
subplot(2,1,2)
semilogx(fv,Make2DFilt(fv, f, 'highpass'),'b')
set(gca,'box','off','tickdir','out','ylim',[0 1.01])
xlabel('Radial frequency (m^{-1})')
ylabel('Filter')
% Filter the DEM
Zhp = SpecFilt2D(Z,dx,dy,f,'highpass');
% Low-pass filter
Zlp = Z - Zhp;
%Zlp = SpecFilt2D(Z,dx,dy,[floLP fhiLP],'lowpass');


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
%Zhp(chan) = mean(Zhp(~chan)); %beautiful!
% Display a shaded relief map of the result (Fig. 7)
%figure(hp)
%ShadePlot(dim.x,dim.y,Zhp)

% Remove edge effects
% You can also see some edge effects produced by the removal of long 
% wavelengths. Let's clip these off.
c = 30;
xclip = c:Nx-c; yclip = c:Ny-c;
Zhp = Zhp(yclip,xclip);
Zlp = Zlp(yclip,xclip);
xhp = dim.x(xclip);
yhp = dim.y(yclip);


% The high-pass-filtered topography
WriteArcGrid(xhp,yhp,Zhp,'orthowv03_ndwi20160707151203_for_wq7_21_highpass_1000_50')
WriteArcGrid(xhp,yhp,Zlp,'orthowv03_ndwi20160707151203_for_wq7_21_highpass_left')
% WriteArcGrid(xhp,yhp,chan,'orthowv01_20150718160144_for_wq7_207_chan')
% WriteArcGrid(xhp,yhp,C,'orthowv01_20150718160144_for_wq7_207_c')