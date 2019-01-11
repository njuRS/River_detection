function outputFileName = spectralanalysis(imageFile,frequency,filterType,inverse,outputpath)

[path,name,ext]=fileparts(imageFile);
fprintf('process %s\n',name);

%kyang, use geotiff files to replace ArcGrid
info = geotiffinfo(imageFile);
[Z, R] = geotiffread(imageFile);

minValue = 0;  % 0 is a safe value for all different data types
Z(isnan(Z)==1)=minValue; %sometimes there are NAN values near image boundaries
Z(isinf(Z)==1)=minValue; 
Z(Z<-1.0)=minValue; %for ndwi image
Z(Z>10000)=minValue; %65535 is the NAN for 16 bit WV images

%Z = single(Z);
Z = single(mat2gray(Z));
dx = info.PixelScale(1);
dy = info.PixelScale(2);
[Ny Nx] = size(Z); % grid dimensions
%for geotiff file write, we cannot change image size, so add new col or row
%and then delete it
signY = 1;
signX = 1;
if(mod(Ny,2)==1)   %fft2D cannot process odd number, so preprocess is required here; debug: 5/19/2016
    newRow = Z(Ny,:);
    Z = [Z ; newRow];
    Ny=Ny + 1;
    signY = 0;
end
if(mod(Nx,2)==1)
    newCol = Z(:,Nx);
    Z = [Z newCol];
    Nx=Nx + 1;
    signX = 0;
end




Z = Detrend(Z);

Zhp = SpecFilt2D(Z,dx,dy,frequency,filterType);


%for geotiff file write, we cannot change image size, so add new col or row
%and then delete it
if(signY==0)   
    Zhp = Zhp(1:Ny-1,:);
end
if(signX==0)
    Zhp = Zhp(:,1:Nx-1);
end

% clip
% c = 1;
% xclip = c:(Nx-c+1); yclip = c:(Ny-c+1);
% Zhp = Zhp(yclip,xclip);
% xhp = dim.x(xclip);
% yhp = dim.y(yclip);

Zhp=mat2gray(Zhp);
Zhp=uint8(Zhp*255);

%invert dark rivers to bright rivers; if the input image is ndwi (in which
%rivers are bright features), use 0
if inverse == 1
    Zhp = 255 - Zhp;
end

outputFileName = [outputpath '\' name '_' filterType ext];
geotiffwrite(outputFileName,Zhp,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
end