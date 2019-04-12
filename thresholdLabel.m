function [rp_img, n_img, fullpath, result_bwlabel] = thresholdLabel(img_num, show_img, method)
% function takes in img_num from single image or stack that you will
% select and loads it, a show_img boolean to define whether you would like to see the
% image at each of 3 steps (load, threshold, bwlabel), and method, the
% thresholding method: Otsu's (1) or based on mean image intensity (2) as
% required by getThreshold. Output is region props of image, number of
% objects in image, fullpath, and the bwlabeled image.

[result, fullpath] = loadAndShowImage(img_num,show_img);
result_1 = getThreshold(result,method);
if show_img == true;
figure, imshow(result_1,[],'InitialMagnification','fit');
end
result_bwlabel = bwlabel(result_1);
if show_img == true;
figure, imshow(result_bwlabel,[],'InitialMagnification','fit');
end
    rp_img = regionprops(result_bwlabel, 'all');
    n_img = length(rp_img);
end
