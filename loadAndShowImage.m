function[result, fullpath]=loadAndShowImage(image_num,show_image,fullpath)
% load and show image, input:(image_num, show_image, optional full_path)=#in stack (integer), true=1/false=0, output opens window to select a .tif stack, loads and shows imagese
if nargin < 3    %check if full_path was given as an argument, if it was not: use uigetfile to select file
[fname, fpath]=uigetfile('*.tif');
fullpath=strcat(fpath,fname);
end
result=imread(fullpath,image_num); %to load image from stack of file that was selected
if (show_image)
    figure
    imshow(result,[],'InitialMagnification','fit');
end  
end
    