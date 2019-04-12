function result = getThreshold(img, method)
% function getThreshold. Input: 2d image (img) and an integer in the range
% 1-2, 1=Otsu's, 2=threshold level is based on mean intensity of image (method)<-- optional, will default to 1. Output: (result) is the thresholded binary or black and white image. 

if nargin < 2
    method = 1
end

if method == 1 %use Otsu's method
    level = graythresh(img);
    result = im2bw(img, level);
elseif method == 2
    img2 = img(:);
    img3 = double(img2);
    level = mean(img2) + 2*std(img3);
    result = img > level;
end

end
