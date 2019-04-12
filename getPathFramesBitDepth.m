
function[fullpath,bitdepth]=getPathFramesBitDepth
% imfinfo to return (1) the full path name of the selected image stack (2) the
%number of spindles and (3) the bit-depth of the first image in an image
%stack whose path you pick interactively, Input: frames_per_spindle,
%Output: fullpath,n_spindles,bitdepth

[fname fpath]=uigetfile('*.tif');
fullpath=strcat(fpath,fname)

info=imfinfo(fullpath); %create a structure with info for each frame
bitdepth=info(1).BitDepth
end




