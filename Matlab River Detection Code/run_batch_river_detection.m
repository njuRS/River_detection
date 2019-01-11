clear;clc;
path = 'E:\Sentinel_2_processing\3_Northwest_GrIS_for_paper_only\ndwi';
outputpath = 'E:\Sentinel_2_processing\3_Northwest_GrIS_for_paper_only\matlab_river_detection';


sensor = 'LandsatNDWI';
%input parameter
inverse = 0;  % set 1 to convert dark rivers into bright rivers (e.g., for panchromatic image); set 0 to keep bright rivers (e.g., for NDWI images)

width = 2;    % small river width for Gabor filter
ppolength = 20; % path opening length
smooth = 0.7;   % smooth paramter for denoise algorithm (denoise algorithm is too slow and scale dependent and thus is abandoned for large images)
histCountThreshold = 100; % a pixel count threhold to stretch image pixel values, default = 1000



cd(path);
files=dir('*.tif');
m=size(files,1);    
for i=1:m
    image = files(i).name;
    if strcmp(sensor,'WV')==1
        f=[1/100 1/20 1/5 1/1];
        filterType='bandpass';
        imageFFT = spectralanalysis(image,f,filterType,inverse,outputpath);
        imageGabor = multidirection_gabor(imageFFT,width,histCountThreshold,outputpath);
        imagePPO = pathopening(imageGabor,ppolength,histCountThreshold,outputpath);   
    elseif strcmp(sensor,'SPOT')==1
        f=[1/200 1/50 1/10 1/5];
        filterType='bandpass';
        imageFFT = spectralanalysis(image,f,filterType,inverse,outputpath);  % inverse = 1, if rivers are dark in the input image
        imageGabor = multidirection_gabor(imageFFT,width,histCountThreshold,outputpath);
        imagePPO = pathopening(imageGabor,ppolength,histCountThreshold,outputpath);
    elseif strcmp(sensor,'SETSM')==1
        f=[1/200 1/100 1/20 1/10];  %follows Karlstrom and Yang, 2016.
        filterType='bandpass';
        smooth = 0.5;
    elseif strcmp(sensor,'Sentinel2')==1
        f=[1/600 1/200 1/40 1/20];
        filterType='bandpass';
        imageFFT = spectralanalysis(image,f,filterType,inverse,outputpath);
        imageGabor = multidirection_gabor(imageFFT,width,histCountThreshold,outputpath);
        imagePPO = pathopening(imageGabor,ppolength,histCountThreshold,outputpath);
    elseif strcmp(sensor,'Landsat')==1
        f=[1/1500 1/600 1/200 1/50];
        filterType='bandpass';
        imageFFT = spectralanalysis(image,f,filterType,inverse,outputpath);
        imageGabor = multidirection_gabor(imageFFT,width,histCountThreshold,outputpath);
        imagePPO = pathopening(imageGabor,ppolength,histCountThreshold,outputpath);
    elseif strcmp(sensor,'LandsatNDWI')==1
        f=[1/1500 1/600 1/200 1/50];
        filterType='bandpass';
        imageFFT = spectralanalysis(image,f,filterType,inverse,outputpath);
        imageGabor = multidirection_gabor(imageFFT,width,histCountThreshold,outputpath);
        imagePPO = pathopening(imageGabor,ppolength,histCountThreshold,outputpath);
    else
        fprintf('do not support this format');
    end
end




