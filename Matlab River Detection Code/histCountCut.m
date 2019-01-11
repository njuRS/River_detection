function output_image = histCountCut( image, threshold )


[hist,edge] = histcounts(image,256);

minBar = 0;
maxBar = 255;

for i=1:256
    if hist(i)<threshold
        continue;
    end
    minBar = i - 1;
    break;
end

for i=1:256
    j = 256 - i;
    if hist(j)<threshold
        continue;
    end
    maxBar = j + 1;
    break;
end

image(image<minBar) = minBar;
image(image>maxBar) = maxBar;

output_image = image;
output_image=mat2gray(output_image);
output_image=uint8(output_image*255);

end

