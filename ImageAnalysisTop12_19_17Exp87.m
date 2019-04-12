clc; 
clear;
close all;
%% SET PARAMETERS
%notes on parametExp87C-3-30miners (check over all parameters before running code!!):
    %for standard analysis of 60x deconvolved MIP of stack - bg (FIJI) use:
        %thresholdMethod = 2;
        %maxarea = 31;
    %for alternative analysis of 60x single frame - bg (FIJI) use:
        %thresholdMethod = 1;
        %maxarea = 101;
    %for analysis of 100x deconvolved MIP of stack - bg (FIJI) CW scope:
        %thresholdMethod = 2;
        %maxarea = 51; %??
        %minarea = 9; %??
   

save = 0; %set to 1 if you want to save figures and summary excel file


DIC = 0; %set to 1 if image for cell counts is DIC and needs to be complemented after loading
clearborder = 0; %set to 1 to clear boarder before thresholding- does not appear to work for Exp 87
optimizeThreshold = 0; %set to 1 if you want to run all thresholding options (NOTE: creates error and does not do downstream analysis), set to 0 to only use threshold method that is defined below
thresholdMethod = 2; % Options: 1, 2, and 3: (1) use arbitrary low threshold of 0.1 (2) Otsu's method (3) greater than mean intensity
fillholes = 0; %set to 1 if using a membrane or cell wall label and hole filling is needed 

% Decide whether to use the watershed processed image (L) or less processed image (I) to ID cells. Set the following variable appropriately

watershed = 0; %set to 1 if you want to use the watershed processed image, 2 if you want to do both

minarea = 2; %not inclusive (cell must be larger than this # of pixels)
maxarea = 50; %not inclusive (cell must be smaller than this # of pixels)--> changed from 31 to 

numOfBins = 20; % for percent histograms

% I = imread(1MATLABanalysis\EfacWGA_DAPI100x\MAX_Exp34-4-007.nd2_decon_ch02.tif);

% \\Client\C$\1MATLABanalysis\EfacWGA_DAPI100x\MAX_Exp34-4-007.nd2_decon_ch02.tif

% I = cell border
% I1 = fluorescence inside cell (DNA or HCR-FISH probes)
% I2 = phase contrast for overlay check

% input number of files, run for loop of the following analysis

files = 1;

%% BIG LOOP
for k = 1:files
    
close all;

%% load and show image, input:(image_num, show_image, optional full_path)=#in stack (integer), true=1/false=0, output opens window to select a .tif stack, loads and shows imagese if nargin < 3    %check if full_path was given as an argument, if it was not: use uigetfile to select file

    [I, fullpath] = loadAndShowImage(1,1);
       title('loaded image');
       
    if save == 1
        fpath = [ fullpath '_analysis\' ];
        mkdir(fpath)
        %saveas(gcf, fullfile(fpath,'loaded image'), 'tif');
        saveas(gcf, fullfile(fpath,'loaded image'), 'fig');
    end
    % [I1, fullpath] = loadAndShowImage(1,1);
    % [I2, fullpath] = loadAndShowImage(1,1);

%% processing for DIC image: complementation... only runs if DIC = 1

    if DIC == 1

    I = imcomplement(I);
   
    figure, imshow(I,[],'InitialMagnification','fit'), title('complemented DIC image');

    end


%% imclearborder to remove cells that touch edge of frame
    if clearborder == 1

    I = imclearborder(I);
    
    figure, imshow(I,[],'InitialMagnification','fit'), title('border cleared');

    end

%% make image black and white
   if thresholdMethod == 1
        I = im2bw(I, 0.1); % Simple method with arbitrary low level for threshold
        figure, imshow(I,[],'InitialMagnification','fit'), title('threshold 0.1');
        if save == 1
            %saveas(gcf, fullfile(fpath,'1 bw arbitrary low thres'), 'tif');
            saveas(gcf, fullfile(fpath,'1 bw arbitrary low thres'), 'fig');
        end
            
   elseif thresholdMethod == 2 
        I = getThreshold(I, 1); %Otsu's method
        figure, imshow(I,[],'InitialMagnification','fit'), title('Otsus method');
        if save == 1
            %saveas(gcf, fullfile(fpath,'2 bw Otsus method'), 'tif');
            saveas(gcf, fullfile(fpath,'2 bw Otsus method'), 'fig');
        end
        
   elseif thresholdMethod == 3 
        I = getThreshold(I, 2); %> mean intensity
        figure, imshow(I,[],'InitialMagnification','fit'), title('>mean intensity');
        if save == 1
            saveas(gcf, fullfile(fpath,'3 bw greater than mean'), 'tif');
            saveas(gcf, fullfile(fpath,'3 bw greater than mean'), 'fig');
        end
   end
    
% Optimization tests    
    if optimizeThreshold == 1

    % Simplest method- no Thresholding (im2bw automatically uses level =
    % 0.1 (arbitrary but very low)
    % where potential range for level is (0,1))... This does not work well,
    % very low % of cells are identified/ pass threshold

    % I_bw = im2bw(I);
    % figure, imshow(I_bw,[],'InitialMagnification','fit');

    % Simple method with arbitrary low level for threshold

    I_bw = im2bw (I, 0.1);
   
    figure, imshow(I_bw,[],'InitialMagnification','fit'), title('threshold 0.1');
    
    if save == 1
    saveas(gcf, fullfile(fpath,'1 bw arbitrary low thres'), 'tif');
    saveas(gcf, fullfile(fpath,'1 bw arbitrary low thres'), 'fig');
    end

    % Otsu's method, holes unfilled
    I_bw1 = getThreshold(I, 1); %Otsu's method
    
    figure, imshow(I_bw1,[],'InitialMagnification','fit'), title('Otsus method');
    
    if save == 1
    saveas(gcf, fullfile(fpath,'2 bw Otsus method'), 'tif');
    saveas(gcf, fullfile(fpath,'2 bw Otsus method'), 'fig');
    end

    % code to bw label objects and display with rainbow colors
    % I_label1 = bwlabel(I_bw1);
    % figure, imshow(I_label1,[],'InitialMagnification','fit');
    % CC = bwconncomp(I_label1);
    % L = labelmatrix(CC);
    % RGB = label2rgb(L);
    % figure, imshow(RGB,[],'InitialMagnification','fit'), title('Otsus method, holes unfilled');


    % % >mean intensity method, holes unfilled
    I_bw2 = getThreshold(I, 2); %> mean intensity
    
    figure, imshow(I_bw2,[],'InitialMagnification','fit'), title('>mean intensity');
    
    if save == 1
    saveas(gcf, fullfile(fpath,'3 bw greater than mean'), 'tif');
    saveas(gcf, fullfile(fpath,'3 bw greater than mean'), 'fig');
    end

    % I_label2 = bwlabel(I_bw2);
    % %figure, imshow(I_label2,[],'InitialMagnification','fit');
    % CC = bwconncomp(I_label2);
    % L = labelmatrix(CC);
    % RGB = label2rgb(L);
    % figure, imshow(RGB,[],'InitialMagnification','fit'), title('>mean intensity, holes unfilled');

    % tests to figure out which im2bw thresholding method to use

    [x,y] = ginput(2); %get clicks for two points from which a line will be drawn
    xline = [x(1), x(2)];
    yline = [y(1), y(1)]; %ignore noise along y axis

    figure
    improfile(I, xline, yline), grid on, title ('arbitrary intensity profile, background subtracted image');
    I_cellinfo = regionprops(I);
    I_cellnumber = length(I_cellinfo);
    if save == 1
    saveas(gcf, fullfile(fpath,'4 arbitrary intensity profile bg sub image'), 'tif');
    saveas(gcf, fullfile(fpath,'4 arbitrary intensity profile bg sub image'), 'fig');
    end

    figure
    improfile(I_bw, xline, yline), grid on, title ('arbitrary intensity profile, arbitrary low thres');
    I_bw_cellinfo = regionprops(I_bw);
    I_bw_cellnumber = length(I_bw_cellinfo);
    if save == 1
    saveas(gcf, fullfile(fpath,'5 arbitrary intensity profile arbitrary low thres'), 'tif');
    saveas(gcf, fullfile(fpath,'5 arbitrary intensity profile arbitrary low thres'), 'fig');
    end

    figure
    improfile(I_bw1, xline, yline), grid on, title ('arbitrary intensity profile, Otsu method');
    I_bw1_cellinfo = regionprops(I_bw1);
    I_bw1_cellnumber = length(I_bw1_cellinfo);
    if save == 1
    saveas(gcf, fullfile(fpath,'6 arbitrary intensity profile Otsus method'), 'tif');
    saveas(gcf, fullfile(fpath,'6 arbitrary intensity profile Otsus method'), 'fig');
    end

    figure
    improfile(I_bw2, xline, yline), grid on, title ('arbitrary intensity profile, >mean intensity');
    I_bw2_cellinfo = regionprops(I_bw2);
    I_bw2_cellnumber = length(I_bw2_cellinfo);
    if save == 1
    saveas(gcf, fullfile(fpath,'7 arbitrary intensity profile greater than mean'), 'tif');
    saveas(gcf, fullfile(fpath,'7 arbitrary intensity profile greater than mean'), 'fig');

    end
    end 
%% fill holes with imfill

    if fillholes == 1

    I = imfill(I, 'holes');
    figure, imshow(I,[],'InitialMagnification','fit'), title('holes filled');

    if save == 1
    saveas(gcf, fullfile(fpath,'2 holes filled'), 'tif');
    saveas(gcf, fullfile(fpath,'2 holes filled'), 'fig');
    end

    I = bwlabel(I);
    %figure, imshow(I_label1A,[],'InitialMagnification','fit');
    CC = bwconncomp(I);
    L = labelmatrix(CC);
    RGB = label2rgb(L);
    figure, imshow(RGB,[],'InitialMagnification','fit'), title('holes filled');

    end

    % probably don't need to "fill holes" unless using a membrane label

%% watershed
     
    % Compute the distance transform of the complement of the binary image.

    D = bwdist(~I);
    figure, imshow(D,[],'InitialMagnification','fit'), title('Distance transform of ~bw');

    % Complement the distance transform, and force pixels that don't belong to the objects to be at -Inf .

    D = -D;
    D(~I) = -Inf;

    % Compute the watershed transform and display the resulting label matrix as an RGB image.

   %L = watershed(D);
    %rgb = label2rgb(L,'winter','w','shuffle');
    %figure, imshow(rgb,'InitialMagnification','fit'), title('Watershed transform of D');


    % L = watershed(D);
    % rgb = label2rgb(L,'jet',[.5 .5 .5]);
    % figure, imshow(rgb,'InitialMagnification','fit'), title('Watershed transform of D')



    %% regionprops to create structure with info about all identified cells

    if watershed == 0 %continue with image that was not processed via watershed
        I_cells = regionprops(I,'all');
        number_cells = length(I_cells);
    elseif watershed == 1 %continue with watershed processed image
        I_cells = regionprops(L,'all');
        number_cells = length(I_cells);
    else %do both
        I_cells = regionprops(I,'all');
        number_cells = length(I_cells);

        I_cells_watershed = regionprops(L,'all');
        number_cells_watershed = length(I_cells);
    end


    %% label cell centroids and plot cell size distribution before filtering objects that are too big to be cells

    % labelc cell centroids
    if watershed == 0
    centroids = cat(1, I_cells(1:end).Centroid);

    figure, imshow(I,'InitialMagnification','fit'), title('I with centroids without watershed');
    hold on
    plot(centroids(:,1),centroids(:,2),'r.');
    hold off
    
    if save == 1
        %saveas(gcf, fullfile(fpath,'I with centroids without watershed'), 'tif');
         saveas(gcf, fullfile(fpath,'I with centroids without watershed'), 'fig');
    end

    elseif watershed == 1
    centroids_watershed = cat(1, I_cells_watershed(2:end).Centroid);

    figure, imshow(I,'InitialMagnification','fit'), title('Watershed transform of D with centroids');
    hold on
    plot(centroids_watershed(:,1),centroids_watershed(:,2),'r.');
    hold off
    
    if save == 1
        saveas(gcf, fullfile(fpath,'Watershed transform of D with centroids'), 'tif');
         saveas(gcf, fullfile(fpath,'Watershed transform of D with centroids'), 'fig');
    end

    else
    centroids = cat(1, I_cells(1:end).Centroid);

    figure, imshow(I,'InitialMagnification','fit'), title('I with centroids without watershed');
    hold on
    plot(centroids(:,1),centroids(:,2),'r.');
    hold off
    
    if save == 1
        saveas(gcf, fullfile(fpath,'I with centroids without watershed'), 'tif');
        saveas(gcf, fullfile(fpath,'I with centroids without watershed'), 'fig');
        
    end

    centroids_watershed = cat(1, I_cells_watershed(2:end).Centroid);

    figure, imshow(I,'InitialMagnification','fit'), title('Watershed transform of D with centroids');
    hold on
    plot(centroids_watershed(:,1),centroids_watershed(:,2),'r.');
    hold off
    
    if save == 1
        saveas(gcf, fullfile(fpath,'Watershed transform of D with centroids'), 'tif');
        saveas(gcf, fullfile(fpath,'Watershed transform of D with centroids'), 'fig');
    end

    end

    % create circles around cells with average diameter (could not get to work)
    % 
    % centers = I_cells(2:end).Centroid;
    % majoraxis = cat(1,I_cells.MajorAxisLength);
    % minoraxis = cat(1,I_cells.MinorAxisLength);
    % diameters = mean([I_cells.(2:end).MajorAxisLength I_cells.(2:end).MinorAxisLength],2);
    % radii = diameters/2;
    % % 
    % figure, imshow(rgb,'InitialMagnification','fit'), title('Watershed transform of D with circles');
    % hold on
    % viscircles(centers,radii);
    % hold off


    % mean area of each object

    
    if watershed == 0
    areas = cat(1,I_cells(1:end).Area);
    meanarea = mean(areas);

    figure, hist(areas), title('cell size distribution in pixels before filtering without watershed');
    
    if save == 1
        %saveas(gcf, fullfile(fpath,'cell size distribution in pixels before filtering without watershed'), 'tif');
        saveas(gcf, fullfile(fpath,'cell size distribution in pixels before filtering without watershed'), 'fig');
    end


    elseif watershed == 1
    areas_watershed = cat(1,I_cells_watershed(1:end).Area);
    meanarea_watershed = mean(areas_watershed);

    figure, hist(areas_watershed), title('cell size distribution in pixels before filtering with watershed');

    if save == 1
        %saveas(gcf, fullfile(fpath,'cell size distribution in pixels before filtering with watershed'), 'tif');
        saveas(gcf, fullfile(fpath,'cell size distribution in pixels before filtering with watershed'), 'fig');
    end
    
    else
    areas = cat(1,I_cells(1:end).Area);
    meanarea = mean(areas);

    areas_watershed = cat(1,I_cells_watershed(1:end).Area);
    meanarea_watershed = mean(areas_watershed);

    figure, hist(areas), title('cell size distribution in pixels before filtering without watershed');
    xlabel('cell size in pixels');
    ylabel('number of cells');
    if save == 1
        %saveas(gcf, fullfile(fpath,'cell size distribution in pixels before filtering without watershed'), 'tif');
        saveas(gcf, fullfile(fpath,'cell size distribution in pixels before filtering without watershed'), 'fig');
    end
    
    figure, hist(areas_watershed), title('cell size distribution in pixels before filtering with watershed');
    xlabel('cell size in pixels');
    ylabel('number of cells');
    if save == 1
        %saveas(gcf, fullfile(fpath,'cell size distribution in pixels before filtering with watershed'), 'tif');
        saveas(gcf, fullfile(fpath,'cell size distribution in pixels before filtering with watershed'), 'fig');
    end
    
    end

    %% eliminate objects that are too big to be cells (greater than 2.5xs mean area)--> (zero centroid and area)
    if watershed == 0
        
        I_label = bwlabel(I); %create bwlabel version of image for filtering purposes
        
        allowableAreaIndexes = (areas > minarea) & (areas < maxarea); %set the minarea and maxarea in first part of code
        keeperIndexes = find(allowableAreaIndexes); %creates list of objects that are in the correct range for cells
        filteredI = ismember(I_label, keeperIndexes); %compare keeper index to objects in unfiltered image, only keep the objects that are in the keeperIndexes
        
        I_cellsfiltered = regionprops(filteredI,'all'); %create structure to describe cells that made it past filter
        number_cells_filtered = length(I_cellsfiltered);
        areasfiltered = cat(1,I_cellsfiltered(1:end).Area);
        meanareafiltered = mean(areasfiltered);
    elseif watershed == 1
        
        L_label = bwlabel(L);
        
        allowableAreaIndexes = (areas > minarea) & (areas < maxarea); %set the minarea and maxarea in first part of code
        keeperIndexes_L = find(allowableAreaIndexes); %creates list of objects that are in the correct range for cells
        filteredL = ismember(L_label, keeperIndexes_L); %compare keeper index to objects in unfiltered image, only keep the objects that are in the keeperIndexes
        
        I_cells_watershedfiltered = regionprops(filteredL,'all'); %create structure to describe cells that made it past filter
        number_cells_watershed_filtered = length(I_cells_watershedfiltered);
        areaswatershedfiltered = cat(1,I_cells_watershedfiltered(1:end).Area);
        meanareawatershedfiltered = mean(areaswatershedfiltered);
        
    else %watershed == 2
        
        %filtering of non-watershed image
        I_label = bwlabel(I); %create bwlabel version of image for filtering purposes
        
        allowableAreaIndexes = (areas > minarea) & (areas < maxarea); %set the minarea and maxarea in first part of code
        keeperIndexes = find(allowableAreaIndexes); %creates list of objects that are in the correct range for cells
        filteredI = ismember(I_label, keeperIndexes); %compare keeper index to objects in unfiltered image, only keep the objects that are in the keeperIndexes
        
        I_cellsfiltered = regionprops(filteredI,'all'); %create structure to describe cells that made it past filter
        number_cells_filtered = length(I_cellsfiltered);
        areasfiltered = cat(1,I_cellsfiltered(1:end).Area);
        meanareafiltered = mean(areasfiltered);
        
        %filtering of watershed image
        L_label = bwlabel(L);
        
        allowableAreaIndexes = (areas > minarea) & (areas < maxarea); %set the minarea and maxarea in first part of code
        keeperIndexes_L = find(allowableAreaIndexes); %creates list of objects that are in the correct range for cells
        filteredL = ismember(L_label, keeperIndexes_L); %compare keeper index to objects in unfiltered image, only keep the objects that are in the keeperIndexes
        
        I_cells_watershedfiltered = regionprops(filteredL,'all'); %create structure to describe cells that made it past filter
        number_cells_watershed_filtered = length(I_cells_watershedfiltered);
        areaswatershedfiltered = cat(1,I_cells_watershedfiltered(1:end).Area);
        meanareawatershedfiltered = mean(areaswatershedfiltered);
    end

    %% label cell centroids and plot cell size distribution after filtering objects that are too big to be cells

    if watershed == 0

        filteredareas = cat(1,I_cellsfiltered(1:end).Area);
        meanfilteredarea = mean(filteredareas);


        filteredcentroids =  cat(1,I_cellsfiltered(1:end).Centroid);

        figure, imshow(I,'InitialMagnification','fit'), title('I with centroids without watershed, filtered');
        xlabel('cell size in pixels');
        ylabel('number of cells');
        hold on
        plot(filteredcentroids(:,1),filteredcentroids(:,2),'r.');
        hold off
        
        if save == 1
        %saveas(gcf, fullfile(fpath,'I with centroids without watershed, filtered'), 'tif');
         saveas(gcf, fullfile(fpath,'I with centroids without watershed, filtered'), 'fig');
        end

        % cell size distribution

        figure, hist(filteredareas), title('cell size distribution without watershed, filtered');
        
        if save == 1
        %saveas(gcf, fullfile(fpath,'cell size distribution without watershed, filtered'), 'tif');
        saveas(gcf, fullfile(fpath,'cell size distribution without watershed, filtered'), 'fig');
        end

        
    elseif watershed == 1

        filteredareas_watershed = cat(1,I_cells_watershedfiltered(1:end).Area);
        meanfilteredarea_watershed = mean(filteredareas_watershed);


        filteredcentroids_watershed =  cat(1,I_cells_watershedfiltered(1:end).Centroid);

        figure, imshow(rgb,'InitialMagnification','fit'), title('Watershed transform of D with centroids, filtered');
        hold on
        plot(filteredcentroids_watershed(:,1),filteredcentroids_watershed(:,2),'r.');
        hold off
        
        if save == 1
        %saveas(gcf, fullfile(fpath,'Watershed transform of D with centroids, filtered'), 'tif');
        saveas(gcf, fullfile(fpath,'Watershed transform of D with centroids, filtered'), 'fig');
        end

        % cell size distribution

        figure, hist(filteredareas_watershed), title('cell size distribution with watershed, filtered');
        xlabel('cell size in pixels');
        ylabel('number of cells');
        
        if save == 1
        %saveas(gcf, fullfile(fpath,'cell size distribution with watershed, filtered'), 'tif');
        saveas(gcf, fullfile(fpath,'cell size distribution with watershed, filtered'), 'fig');
        end

    else

        filteredareas = cat(1,I_cellsfiltered(1:end).Area);
        meanfilteredarea = mean(filteredareas);
        
        filteredcentroids =  cat(1,I_cellsfiltered(1:end).Centroid);

        figure, imshow(rgb,'InitialMagnification','fit'), title('I with centroids without watershed, filtered');
        hold on
        plot(filteredcentroids(:,1),filteredcentroids(:,2),'r.');
        hold off
        
        if save == 1
        %saveas(gcf, fullfile(fpath,'I with centroids without watershed, filtered'), 'tif');
        saveas(gcf, fullfile(fpath,'I with centroids without watershed, filtered'), 'fig');
        end

        % cell size distribution

        figure, hist(filteredareas), title('cell size distribution without watershed, filtered');
        xlabel('cell size in pixels');
        ylabel('number of cells');
        
        if save == 1
        %saveas(gcf, fullfile(fpath,'cell size distribution without watershed, filtered'), 'tif');
        saveas(gcf, fullfile(fpath,'cell size distribution without watershed, filtered'), 'fig');
        end

        filteredareas_watershed = cat(1,I_cells_watershedfiltered(1:end).Area);
        meanfilteredarea_watershed = mean(filteredareas_watershed);


        filteredcentroids_watershed =  cat(1,I_cells_watershedfiltered(1:end).Centroid);

        figure, imshow(rgb,'InitialMagnification','fit'), title('Watershed transform of D with centroids, filtered');
        hold on
        plot(filteredcentroids_watershed(:,1),filteredcentroids_watershed(:,2),'r.');
        hold off
        
        if save == 1
        %saveas(gcf, fullfile(fpath,'Watershed transform of D with centroids, filtered'), 'tif');
        saveas(gcf, fullfile(fpath,'Watershed transform of D with centroids, filtered'), 'fig');
        end

        % cell size distribution

        figure, hist(filteredareas_watershed), title('cell size distribution with watershed, filtered');
        
        if save == 1
       % saveas(gcf, fullfile(fpath,'cell size distribution with watershed, filtered'), 'tif');
        saveas(gcf, fullfile(fpath,'cell size distribution with watershed, filtered'), 'fig');
        end

    end

            %% overlay filtered and unfiltered cells on original image

            figure, imshow(I,'InitialMagnification','fit'), title('Watershed transform of D with centroids, original image');
            hold on
            plot(centroids(:,1),centroids(:,2),'ro');
            plot(filteredcentroids(:,1),filteredcentroids(:,2),'y.');
            hold off
            
            if save == 1
         %       saveas(gcf, fullfile(fpath,'Watershed transform of D with centroids, original image'), 'tif');
                saveas(gcf, fullfile(fpath,'Watershed transform of D with centroids, original image'), 'fig');
            end

            %% overlay other cell id options

%             I_label1rp = regionprops(I_label1);
% 
%             centroids_label1rp = cat(1, I_label1rp(2:end).Centroid);
% 
%             hold on
%             plot(centroids_label1rp(:,1),centroids_label1rp(:,2),'co'), plot(centroids(:,1),centroids(:,2),'ro');
%             plot(filteredcentroids(:,1),filteredcentroids(:,2),'y.');
%             hold off
%             
%             if save == 1
%                 saveas(gcf, fullfile(fpath,'Watershed transform of D with centroids, original image, more'), 'tif');
%             end


            %% imoverlay function: doesn't seem to be working right
            % edges = edge(L);
            % new1 = imoverlay(I,edges,[1 0 1]);
            % figure, imshow(new1,'InitialMagnification','fit')


            %% load fluorescence image to add intensity data to current structure

            [I1, fullpath] = loadAndShowImage(1,1);

    %% Add fluoresecence field (HCR-FISH) to identified cells structure 
        % unfiltered data
    for i = 1:number_cells

        cellpixellocations = I_cells(i).PixelList(:,:);
        dim = size(cellpixellocations);
        I_cells(i).I1intensity = zeros(dim(1),1);

        for j = 1:dim(1)

            I_cells(i).I1intensity(j) = I1(cellpixellocations(j,2),cellpixellocations(j,1));
        end

    end    
        
    %% Find average fluorescence from HCR-FISH channel for each cell and add to identified cells structure

    % unfiltered data
    for i = 1:number_cells

        I_cells(i).I1averageintensity = mean(I_cells(i).I1intensity);

    end

    I1averageintensity =  cat(1,I_cells(1:end).I1averageintensity);

    % for i = 1:number_cells
    %     
    %     if I1averageintensity(i) >100 
    %         
    %         I1averageintensity(i) = NaN;
    %         
    %     end
    %     
    % end

    figure, hist(I1averageintensity), title('histogram of average fluorescent intensity in each cell');
    xlabel('mean fluorescent intensity (AU)');
    ylabel('number of cells');
    axis([0,(8*10^2),0,500]);
    
    if save == 1
        % saveas(gcf, fullfile(fpath,'histogram of average HCR-FISH fluorescent intensity in identified cells (without watershed)'), 'tif');
         saveas(gcf, fullfile(fpath,'histogram of average HCR-FISH fluorescent intensity in identified cells (without watershed)'), 'fig');
    end
    
    
    [histFreq, histXout] = hist(I1averageintensity, numOfBins);
    figure;
    bar(histXout, histFreq/sum(histFreq)*100), title('percent of cells in average fluorescent intensity bins');
    xlabel('mean fluorescent intensity (AU)');
    ylabel('percent of cells');
    axis([0,(8*10^2),0,100]);

    if save == 1
       % saveas(gcf, fullfile(fpath, 'percent of cells in average fluorescent intensity bins'), 'tif');
        saveas(gcf, fullfile(fpath, 'percent of cells in average fluorescent intensity bins'), 'fig');
    end
    
    
    %% Find sum of fluorescence from HCR-FISH channel for each cell and add to identified cells structure

    % unfiltered data
    for i = 1:number_cells

        I_cells(i).I1sumintensity = sum(I_cells(i).I1intensity);
    end

    I1sumintensity = cat(1,I_cells(1:end).I1sumintensity);

    figure, hist(I1sumintensity), title('histogram of sum fluorescent intensity in each cell');
    xlabel('sum fluorescent intensity in each cell (AU)');
    ylabel('number of cells');
    axis([0,(5*10^3),0,500]);
    
    if save == 1
        % saveas(gcf, fullfile(fpath,'histogram of sum HCR-FISH fluorescent intensity in identified cells (without watershed)'), 'tif');
         saveas(gcf, fullfile(fpath,'histogram of sum HCR-FISH fluorescent intensity in identified cells (without watershed)'), 'fig');
    end


    %% plot cells intensity (with standard deviation and without)
    
        % unfiltered data
    for i = 1:number_cells

        I_cells(i).standardDev = std(I_cells(i).I1intensity);

    end

    standardDev = cat(1,I_cells(1:end).standardDev);


    % figure, scatter(I1averageintensity, number_cells);
    figure, errorbar(I1averageintensity, standardDev, 'rx'), title('average fluorescence and st dev per cell');
    xlabel('cell number');
    ylabel('mean fluorescent intensity (AU) with standard deviation');
    axis([0,400,(-5000),(8*10^2)]);
    
    if save == 1
         %saveas(gcf, fullfile(fpath,'average HCR-FISH fluorescent intensity and standard dev'), 'tif');
         saveas(gcf, fullfile(fpath,'average HCR-FISH fluorescent intensity and standard dev'), 'fig');
    end
    
    figure, scatter((1:number_cells),I1averageintensity, 'rx'), title('average fluorescence per cell'); %added on 5/12/2016
    xlabel('cell number');
    ylabel('mean fluorescent intensity (AU)');
    axis([0,400,0,(8*10^2)]);
    
    if save == 1
        %saveas(gcf, fullfile(fpath, 'average HCR-FISH fluorescent intensity'), 'tif'); 
        saveas(gcf, fullfile(fpath, 'average HCR-FISH fluorescent intensity'), 'fig'); 
    end

    figure, errorbar(I1sumintensity, standardDev, 'rx'), title('sum fluorescence and st dev');
    xlabel('cell number');
    ylabel('sum fluorescent intensity (AU) with standard deviation of mean');
    axis([0,400,0,(5*10^3)]);
    
    if save == 1
       %  saveas(gcf, fullfile(fpath,'sum HCR-FISH fluorescent intensity and standard dev'), 'tif');
         saveas(gcf, fullfile(fpath,'sum HCR-FISH fluorescent intensity and standard dev'), 'fig');
    end
    
 
    %% plot cell size (area) vs. sum intensity (HCR-FISH)

    %not working yet
    %figure, plot(areas, I1sumintensity(2:end)), title('sum HCR-FISH fluorescent intensity vs. cell size (area)');


    % rp = round([rp_dapi(i).Centroid]); %determine the center of each dapi object
    %     if result_agg2(rp(1,2), rp(1,1)) ~= 0
    %         dapi_agg_overlap_data(i,1) = 1; %dapi puncta overlaps with agg at this mean centroid location
    %     else
    %         dapi_agg_overlap_data(i) = 0; %no overlap
    %     end
    %     dapi_agg_overlap_data(i,2) = rp(1,2); %pixel row location
    %     dapi_agg_overlap_data(i,3) = rp(1,1); %pixel column location
    %     dapi_agg_overlap_data(i,4) = result_agg2(rp(1,2),rp(1,1)); % w
    % 

    
    
    %% Add fluoresecence field (HCR-FISH) to identified cells structure (filtered)

    for i = 1:number_cells_filtered

        cellpixellocations = I_cellsfiltered(i).PixelList(:,:);
        dim = size(cellpixellocations);
        I_cellsfiltered(i).I1intensity = zeros(dim(1),1);

        for j = 1:dim(1)

            I_cellsfiltered(i).I1intensity(j) = I1(cellpixellocations(j,2),cellpixellocations(j,1));
        end

    end

    %% Find average fluorescence from HCR-FISH channel for each cell and add to identified cells structure (filtered)

    for i = 1:number_cells_filtered

        I_cellsfiltered(i).I1averageintensity = mean(I_cellsfiltered(i).I1intensity);

    end

    I1averageintensityfiltered =  cat(1,I_cellsfiltered(1:end).I1averageintensity);

    % for i = 1:number_cells
    %     
    %     if I1averageintensity(i) >100 
    %         
    %         I1averageintensity(i) = NaN;
    %         
    %     end
    %     
    % end

    figure, hist(I1averageintensityfiltered), title('histogram of average fluorescent intensity in each cell- filtered');
    xlabel('mean fluorescent intensity (AU)');
    ylabel('number of cells');
    axis([0,(8*10^2),0,400]);
    
    if save == 1
      %   saveas(gcf, fullfile(fpath,'histogram of average HCR-FISH fluorescent intensity in identified cells (without watershed) filtered'), 'tif');
         saveas(gcf, fullfile(fpath,'histogram of average HCR-FISH fluorescent intensity in identified cells (without watershed) filtered'), 'fig');
    end
    
    [histFreq, histXout] = hist(I1averageintensityfiltered, numOfBins);
    figure;
    bar(histXout, histFreq/sum(histFreq)*100), title('percent of cells in average fluorescent intensity bins- filtered');
    xlabel('mean fluorescent intensity (AU)');
    ylabel('percent of cells');
    axis([0,(8*10^2),0,100]);

    if save == 1
     %   saveas(gcf, fullfile(fpath, 'percent of cells in average fluorescent intensity bins- filtered'), 'tif');
        saveas(gcf, fullfile(fpath, 'percent of cells in average fluorescent intensity bins- filtered'), 'fig');
    end

    %% Find sum of fluorescence from HCR-FISH channel for each cell and add to identified cells structure (filtered)

    for i = 1:number_cells_filtered

        I_cellsfiltered(i).I1sumintensity = sum(I_cellsfiltered(i).I1intensity);
    end

    I1sumintensityfiltered = cat(1,I_cellsfiltered(1:end).I1sumintensity);

    figure, hist(I1sumintensityfiltered), title('histogram of sum fluorescent intensity in each cell- filtered');
    xlabel('sum fluorescent intensity in each cell (AU)');
    ylabel('number of cells');
    axis([0,(5*10^3),0,400]);
    
    if save == 1
    %     saveas(gcf, fullfile(fpath,'histogram of sum HCR-FISH fluorescent intensity in identified cells (without watershed) filtered'), 'tif');
         saveas(gcf, fullfile(fpath,'histogram of sum HCR-FISH fluorescent intensity in identified cells (without watershed) filtered'), 'fig');
    end

    %%
    for i = 1:number_cells_filtered

        I_cellsfiltered(i).standardDev = std(I_cellsfiltered(i).I1intensity);

    end

    standardDevfiltered = cat(1,I_cellsfiltered(1:end).standardDev);


    figure, errorbar(I1averageintensityfiltered, standardDevfiltered, 'rx'), title('average fluorescence and st dev per cell- filtered');
    xlabel('cell number');
    ylabel('mean fluorescent intensity (AU) with standard deviation');
    axis([0,400,(-5000),(8*10^2)]);
    
    if save == 1

   % saveas(gcf, fullfile(fpath,'average HCR-FISH fluorescent intensity and standard dev filtered'), 'tif');
    saveas(gcf, fullfile(fpath,'average HCR-FISH fluorescent intensity and standard dev filtered'), 'fig');
    
    end
    
    figure, scatter((1:number_cells_filtered),I1averageintensityfiltered, 'rx'), title('average HCR-FISH fluorescent intensity in each cell'); %added on 5/12/2016
    xlabel('cell number');
    ylabel('mean fluorescent intensity (AU)');
    axis([0,400,0,(8*10^2)]);
    
    if save == 1
        %saveas(gcf, fullfile(fpath, 'average HCR-FISH fluorescent intensity filtered'), 'tif'); 
        saveas(gcf, fullfile(fpath, 'average HCR-FISH fluorescent intensity filtered'), 'fig'); 
    end
    

    figure, errorbar(I1sumintensityfiltered, standardDevfiltered, 'rx'), title('sum HCR-FISH fluorescent intensity and standard dev filtered');
    xlabel('cell number');
    ylabel('sum fluorescent intensity (AU) with standard deviation of mean');
    axis([0,400,0,(5*10^3)]);
    
    if save == 1
         %saveas(gcf, fullfile(fpath,'sum HCR-FISH fluorescent intensity and standard dev filtered'), 'tif');
         saveas(gcf, fullfile(fpath,'sum HCR-FISH fluorescent intensity and standard dev filtered'), 'fig');
    end
    

    %% plot cell size (area) vs. sum intensity (HCR-FISH)

    %not working yet
    %figure, plot(areas, I1sumintensity(2:end)), title('sum HCR-FISH fluorescent intensity vs. cell size (area)');


    % rp = round([rp_dapi(i).Centroid]); %determine the center of each dapi object
    %     if result_agg2(rp(1,2), rp(1,1)) ~= 0
    %         dapi_agg_overlap_data(i,1) = 1; %dapi puncta overlaps with agg at this mean centroid location
    %     else
    %         dapi_agg_overlap_data(i) = 0; %no overlap
    %     end
    %     dapi_agg_overlap_data(i,2) = rp(1,2); %pixel row location
    %     dapi_agg_overlap_data(i,3) = rp(1,1); %pixel column location
    %     dapi_agg_overlap_data(i,4) = result_agg2(rp(1,2),rp(1,1)); % w
    % 
    %% compile data from this field and other fields
    
    if files > 1
    
    %structure
    if k == 1
        I_cells_fields = I_cellsfiltered; %initiate compiled structure in i=1
        
    else 
        I_cells_fields_temp = [I_cellsfiltered; I_cells_fields]; %combine cell info structure from current field with previous fields
        I_cells_fields = I_cells_fields_temp; %carry compiled data to next i in loop
    end
    
    %cell size
    meanareafiltered;
    if k == 1
        meanareas = meanareafiltered; %initiate compiled table in i=1
    else 
        meanareas_temp = [meanareafiltered; meanareas]; %combine cell info from current field with previous fields
        meanareas = meanareas_temp; 
    end
    
    %mean intensity
    averagemeanintensity = mean(I1averageintensityfiltered);
    if k == 1 
        averagemeanintensities = averagemeanintensity; 
    else 
        averagemeanintensity_temp = [averagemeanintensity; averagemeanintensities];
        averagemeanintensities = averagemeanintensity_temp;
    end
    
    %sum intensity
    averagesumintensity = mean(I1sumintensityfiltered);
    if k == 1
        averagesumintensities = averagesumintensity;
    else 
        averagesumintensity_temp = [averagesumintensity; averagesumintensities];
        averagesumintensities = averagesumintensity_temp;
    end
    
    %standard dev
    averagestandarddev = mean(standardDevfiltered);
    if k == 1
        averagestandarddevs = averagestandarddev;
    else
        averagestandarddev_temp = [averagestandarddev; averagestandarddevs];
        averagestandarddevs = averagestandarddev_temp;
    end
    
    end
    
%% compile summary data to save
majoraxislengths = cat(1,I_cellsfiltered.MajorAxisLength);
minoraxislengths = cat(1,I_cellsfiltered.MinorAxisLength);

    summary = zeros(number_cells_filtered(end),10); %initiate matrix for compiled cell data to export, column 1 for field ID
    if files > 1
        %fieldnumber(1:number_cells_filtered,1) = k; % I don't think I need
        %this if the following line "summary(1:end,1) = fieldnumber;" has
        %fieldnumber switched to k. ie: "summary(1:end,1) = k;"
        summary(1:end,1) = k;
    end
    summary(1:end,2) = 1:number_cells_filtered(end); %cell number
    summary(1:end,3) = I1averageintensityfiltered; %average intensity
    summary(1:end,4) = I1sumintensityfiltered; %sum intensity
    summary(1:end,5) = standardDevfiltered; %standard deviation of intensity in each cell
    summary(1:end,6) = filteredareas; %number of pixels in each cell
    summary(1:end,7:8) = filteredcentroids; %coordinates of the centroid of each cell
    summary(1:end,9) = majoraxislengths; 
    summary(1:end,10) = minoraxislengths;
    
if files > 1

    if k == 1
    fieldssummary = summary;
    fieldssummary(1:number_cells_filtered,1) = 1;
    
    elseif k > 1
        
        fieldssummary = [fieldssummary; summary]; %need to concatenate fieldssummary with summary (this will add the most recent field to the previous fields summary)
        
    end 
end    
end %end of big 'for loop'

%% save data- outside big loop
if save == 1
if files == 1
    xlswrite(fullfile(fpath,'summary'),summary); %save summary
elseif files > 1
    xlswrite(fullfile(fpath, 'fieldssummary'), fieldssummary); %save fieldssummary 
end
end
%% plot data from multiple fields

if files > 1
% plot data from cells from multiple fields

fieldsaverageintensity =  cat(1,I_cells_fields(1:end).I1averageintensity);

    figure, hist(fieldsaverageintensity), title('histogram of average fluorescent intensity in each cell- multiple filtered fields');
    xlabel('mean fluorescent intensity (AU)');
    ylabel('number of cells');
    axis([0,(8*10^4), 0, inf]); %set x-axis but not y-axis... it will default to max # of cells in that bin
    %axis([0,(5*10^4),0,500]);
    
    if save == 1
         %saveas(gcf, fullfile(fpath,'1histogram of average HCR-FISH fluorescent intensity in identified cells- multiple fields (without watershed) filtered'), 'tif');
         saveas(gcf, fullfile(fpath,'1histogram of average HCR-FISH fluorescent intensity in identified cells- multiple fields (without watershed) filtered'), 'fig');
    end

    [histFreq, histXout] = hist(fieldsaverageintensity, numOfBins);
    figure;
    bar(histXout, histFreq/sum(histFreq)*100), title('percent of cells in average fluorescent intensity bins- multiple filtered fields');
    xlabel('mean fluorescent intensity (AU)');
    ylabel('percent of cells');
    axis([0,(8*10^4),0,100]);

    if save == 1
        %saveas(gcf, fullfile(fpath, '1percent of cells in average fluorescent intensity bins- multiple fields filtered'), 'tif');
        saveas(gcf, fullfile(fpath, '1percent of cells in average fluorescent intensity bins- multiple fields filtered'), 'fig');
    end
    
    
fieldssumintensity = cat(1,I_cells_fields(1:end).I1sumintensity);

    figure, hist(fieldssumintensity), title('histogram of sum fluorescent intensity in each cell- multiple filtered fields');
    xlabel('sum fluorescent intensity in each cell (AU)');
    ylabel('number of cells');
    axis([0,(5*10^6),0,inf]);
    %axis([0,(5*10^6),0,200]);
    
    if save == 1
         %saveas(gcf, fullfile(fpath,'1histogram of sum HCR-FISH fluorescent intensity in identified cells- multiple fields (without watershed)filtered'), 'tif');
         saveas(gcf, fullfile(fpath,'1histogram of sum HCR-FISH fluorescent intensity in identified cells- multiple fields (without watershed)filtered'), 'fig');
    end

fieldsstandardDev = cat(1,I_cells_fields(1:end).standardDev);


    figure, errorbar(fieldsaverageintensity, fieldsstandardDev, 'rx'), title('average fluorescence and st dev per cell- multiple filtered fields');
    xlabel('cell number');
    ylabel('mean fluorescent intensity (AU) with standard deviation');
    axis([0,inf,(-5000),(8*10^2)]); 
    %axis([0,500,(-2000),(6*10^4)]);
    
    if save == 1
         %saveas(gcf, fullfile(fpath,'1average HCR-FISH fluorescent intensity and standard dev- multiple fields filtered'), 'tif');
         saveas(gcf, fullfile(fpath,'1average HCR-FISH fluorescent intensity and standard dev- multiple fields filtered'), 'fig');
    end
    
    number_cells_fields = length(fieldssummary(:,1)); 
    figure, scatter((1:number_cells_fields),fieldsaverageintensity, 'rx'), title('average fluorescent intensity in each cell- multiple fields'); %added on 5/12/2016
    xlabel('cell number');
    ylabel('mean fluorescent intensity (AU)');
    axis([0,inf,0,(8*10^2)]);
    
    if save == 1
        %saveas(gcf, fullfile(fpath, '1average HCR-FISH fluorescent intensity multiple filtered fields'), 'tif'); 
        saveas(gcf, fullfile(fpath, '1average HCR-FISH fluorescent intensity multiple filtered fields'), 'fig'); 
    end
    
    figure, errorbar(fieldssumintensity, fieldsstandardDev, 'rx'), title('sum HCR-FISH fluorescent intensity and st dev- filtered fields');
    xlabel('cell number');
    ylabel('sum fluorescent intensity (AU) with standard deviation of mean');
    axis([0,inf,0,(5*10^3)]);
    %axis([0,500,0,(1*10^7)]);
    
    if save == 1
         %saveas(gcf, fullfile(fpath,'1sum HCR-FISH fluorescent intensity and standard dev- multiple fields filtered'), 'tif');
         saveas(gcf, fullfile(fpath,'1sum HCR-FISH fluorescent intensity and standard dev- multiple fields filtered'), 'fig');
    end 
    
    % plot averages from each field to compare field variation
    
    figure, bar(meanareas), title('mean cell area for each field filtered');
    xlabel('cell field');
    ylabel('mean cell area (pixels)');
    
    if save == 1
         %saveas(gcf, fullfile(fpath,'1mean cell area for each field filtered'), 'tif');
         saveas(gcf, fullfile(fpath,'1mean cell area for each field filtered'), 'fig');
    end
    
    figure, bar(averagemeanintensities), title('average mean intensities of cells in each field filtered ');
    xlabel('cell field');
    ylabel('average of mean cell intensities (AU)');
    
    if save == 1
         %saveas(gcf, fullfile(fpath,'1average mean intensities of cells in each field filtered'), 'tif');
         saveas(gcf, fullfile(fpath,'1average mean intensities of cells in each field filtered'), 'fig');
    end
    
    figure, bar(averagesumintensities), title('average sum intensities of cells in each field filtered');
    xlabel('cell field');
    ylabel('average of sum cell intensity (AU)');
    
    if save == 1
         %saveas(gcf, fullfile(fpath,'1average sum intensities of cells in each field filtered'), 'tif');
         saveas(gcf, fullfile(fpath,'1average sum intensities of cells in each field filtered'), 'fig');
    end
    
    figure, bar(averagestandarddevs), title('average standard dev of fluor intensity of cells in each field filtered');
    xlabel('cell field');
    ylabel('average of standard dev of cell intensity');
    
    if save == 1
         %saveas(gcf, fullfile(fpath,'1average standard dev of fluor intensity of cells in each field filtered'), 'tif');
         saveas(gcf, fullfile(fpath,'1average standard dev of fluor intensity of cells in each field filtered'), 'fig');
    end
    
end

    