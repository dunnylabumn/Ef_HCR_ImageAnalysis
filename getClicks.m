function[col, row]= getClicks(my_image,show_clicks)
%getClicks takes an image (my_image=2d or 3d matrix) and an optional
%Boolean (show_clicks), displays image, allows clicks, returns col row
%locations of clicks. 
if nargin < 2
    show_clicks=true;
end

imshow(my_image,[],'InitialMagnification','fit'); % display image
text(10,10,'Click on desired locations then hit enter','Color','white') %take in clicks

[col, row, P]=impixel(); %return the column (col) and row (row) locations of the clicks
close; %close image

 if show_clicks
     imshow(my_image,[],'InitialMagnification','fit');
     hold on;
     scatter(col, row);
     hold off;
 end
end