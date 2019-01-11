function outputFileName = denoise(imageFile,smoothParameter,outputpath)

[path,name,ext]=fileparts(imageFile);
fprintf('process %s\n',name);
%params
M = 7;      % search area size (2*M + 1)^2
alpha = 3;  % patch size (2*alpha + 1)^2
h = smoothParameter;    % smoothing parameter [0-infinite]. Kang, 0.7 for WV & SPOT images; 0.5 for SETSM
% If you can see structures in the "Residual image" decrease this parameter
% If you can see speckle in the "denoised image" increase this parameter
% h = 0.3; % smoothing parameter for Landsat image only 
offset = 100; % to avoid Nan in Pearson divergence computation
% According to the gain used by your US device this offset can be adjusted.
info = geotiffinfo(imageFile);
[img,R] = geotiffread(imageFile);

[dimxy dimt] = size(size(img));
if ( dimt > 2)
     img = rgb2gray(img);
end

% Intensity normalization
imgd = double(img);
mini = (min(imgd(:)));
imgd = (imgd - mini);
maxi = max(imgd(:));
imgd = (imgd / maxi) * 255;
imgd = imgd + offset; % add offset to enable the pearson divergence computation (i.e. avoid division by zero).
s = size(imgd);

% Padding
imgd = padarray(imgd,[alpha alpha],'symmetric');
fimgd=bnlm2D(imgd,M,alpha,h);
fimgd = fimgd - offset;
imgd = imgd - offset;
imgd = imgd(alpha+1: s(1)+alpha, alpha+1: s(2)+alpha);
fimgd = fimgd(alpha+1: s(1)+alpha, alpha+1: s(2)+alpha);

fimg = fimgd/ max(imgd(:));

%fimg=max(fimg(:))-fimg; %invert dark rivers to bright


fimg=mat2gray(fimg);
fimg=uint8(fimg*255);


outputFileName=[outputpath '\' name '_norm_denoised' ext];
geotiffwrite(outputFileName,fimg,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);