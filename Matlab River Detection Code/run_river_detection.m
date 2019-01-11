clear;clc;
path = 'C:\Users\lenovo\Desktop';
image = [path '\' 'test_sentinel2_image.tif'];
outputpath = path;
sensor = 'Sentinel2';
inverse = 0; % set 1 to convert dark rivers into bright rivers (e.g., for panchromatic image); set 0 to keep bright rivers (e.g., for NDWI images)

% path = 'E:\Sentinel_2_processing\Landsat_comparison';
% image = [path '\' 'ndwi_pansharpen_stack_LC80070132016208LGN00.tif'];
% outputpath = 'E:\Sentinel_2_processing\Landsat_comparison';
% sensor = 'Landsat';
% inverse = 0; % set 1 to convert dark rivers into bright rivers (e.g., for panchromatic image); set 0 to keep bright rivers (e.g., for NDWI images)

width = 2;
ppolength = 20;
smooth = 0.7;
histCountThreshold = 1000;

if strcmp(sensor,'WV')==1  % 0.5 m resolution
    f=[1/100 1/20 1/5 1/1];
    filterType='bandpass';
    imageFFT = spectralanalysis(image,f,filterType,inverse,outputpath);
    imageGabor = multidirection_gabor(imageFFT,width,histCountThreshold,outputpath);
    imagePPO = pathopening(imageGabor,ppolength,histCountThreshold,outputpath);   
elseif strcmp(sensor,'SPOT')==1 % 1.5 m resolution
    f=[1/200 1/50 1/10 1/5];
    filterType='bandpass';
    imageFFT = spectralanalysis(image,f,filterType,inverse,outputpath);  % inverse = 1, if rivers are dark in the input image
    imageGabor = multidirection_gabor(imageFFT,width,histCountThreshold,outputpath);
    imagePPO = pathopening(imageGabor,ppolength,histCountThreshold,outputpath);
elseif strcmp(sensor,'SETSM')==1 % 2.0 m resolution
    f=[1/200 1/100 1/20 1/10];  %follows Karlstrom and Yang, 2016.
    filterType='bandpass';
    imageFFT = spectralanalysis(image,f,filterType,inverse,outputpath);  % inverse = 1, if rivers are dark in the input image
    imageGabor = multidirection_gabor(imageFFT,width,histCountThreshold,outputpath);
    imagePPO = pathopening(imageGabor,ppolength,histCountThreshold,outputpath);
elseif strcmp(sensor,'Sentinel2')==1 % 10 m resolution
    f=[1/600 1/200 1/40 1/20];
    filterType='bandpass';
    imageFFT = spectralanalysis(image,f,filterType,inverse,outputpath);
    imageGabor = multidirection_gabor(imageFFT,width,histCountThreshold,outputpath);
    imagePPO = pathopening(imageGabor,ppolength,histCountThreshold,outputpath);
elseif strcmp(sensor,'Landsat')==1 % 15 m resolution
    f=[1/1000 1/500 1/200 1/50];
    filterType='bandpass';
    imageFFT = spectralanalysis(image,f,filterType,inverse,outputpath);
    imageGabor = multidirection_gabor(imageFFT,width,histCountThreshold,outputpath);
    imagePPO = pathopening(imageGabor,ppolength,histCountThreshold,outputpath);
elseif strcmp(sensor,'LandsatNDWI')==1 % 30 m resolution
    f=[1/1000 1/500 1/200 1/50];
    filterType='bandpass';
    imageFFT = spectralanalysis(image,f,filterType,inverse,outputpath);
    imageGabor = multidirection_gabor(imageFFT,width,histCountThreshold,outputpath);
    imagePPO = pathopening(imageGabor,ppolength,histCountThreshold,outputpath);
else
    fprintf('do not support this format');
end
