path = 'E:\0000_processed_Sentinel2\southwest_GrIS\matlab_river_detection\part1';
imageFile = [path '\' 'ndwi_S2A_OPER_MSI_L1C_TL_SGS__20160804T202340_A005842_T22WFA_bandpass_gabor_cpo20.tif'];

[path,name,ext]=fileparts(imageFile);
fprintf('process %s\n',name);
info = geotiffinfo(imageFile);
[image,R] = geotiffread(imageFile);

image = image>3.71;

image = bwlabel(image);

outputpath = 'E:';

outputFileName = [outputpath '\' name '_river' ext];
geotiffwrite(outputFileName,image,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);