function outputFileName = pathopening(imageFile,lengthThreshold,histCountThreshold,outputpath)

[path,name,ext]=fileparts(imageFile);
fprintf('process %s\n',name);

%initialize dip_image toolbox
run('C:\Program Files\DIPimage 2.8.1\dipstart.m');

info = geotiffinfo(imageFile);
[image,R] = geotiffread(imageFile);

image=single(image);

%I planed to use mask image to reduce computational time but it seems this
%mask image does not work in path opening function
mask = ones(size(image));

dip_mask = dip_image(mask,'bin');
dip_data = dip_image(image,'dfloat');

out = dip_pathopening(dip_data,dip_mask,lengthThreshold,0,1);
image_path_opened = dip_array(out);

%pixels smaller than mean cannot be rivers. This stretch is used for final
%display purpose
meanPO = mean(image_path_opened(:));
image_path_opened(image_path_opened<meanPO) = meanPO;

image_path_opened=mat2gray(image_path_opened);
image_path_opened=uint8(image_path_opened*255);

image_path_opened = histCountCut(image_path_opened,histCountThreshold);

outputFileName=[outputpath '\' name '_cpo' num2str(lengthThreshold) ext];
geotiffwrite(outputFileName,image_path_opened,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);

end