function outputFileName = multidirection_gabor(imageFile,width,histCountThreshold,outputpath,varargin)
[path,name,ext]=fileparts(imageFile);
fprintf('process %s\n',name);

info = geotiffinfo(imageFile);
[image,R] = geotiffread(imageFile);
image=single(image);

if nargin == 4
    %default 
    elongation = 1.0; %determines the elongation of the Gabor filter in the orientation direction, with respect to its thickness.
    filterType = 'even';
elseif nargin == 5
    elongation = cell2mat(varargin(1));
    filterType = 'even';
elseif nargin == 6
    elongation = cell2mat(varargin(1));
    filterType = cell2mat(varargin(2));
else
    fprintf('incorret input parameters\n');
end
    

angles = 0:15:165;
[l w] = size(image);
%g_total = zeros(l,w,12);
image_gf = zeros(l,w);  %avoid to create large matrix;

for i = 1:12
    theta = angles(i)*(pi/180); %angle of rotation
    gamma = 2 * width;
    lambda = 2.0 * width; %1/f
    psi = 0; %offset
    sigma = gamma/(2*sqrt(2*log(2)));
    gb_filt=gaborfilter(sigma,theta,lambda,psi,elongation,filterType);
    g = conv2(image,gb_filt,'same');
    image_gf = max(image_gf,g);
    %g_total(:,:,i) = g;
end

%image_gf = max(g_total,[],3);
image_gf=mat2gray(image_gf);
image_gf=uint8(image_gf*255);


image_gf = histCountCut(image_gf,histCountThreshold);

if strcmp(filterType, 'even')
    outputFileName=[outputpath '\' name '_gabor' ext];
else
    outputFileName=[outputpath '\' name '_gabor_odd' ext];
end

geotiffwrite(outputFileName,image_gf,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);

end

    