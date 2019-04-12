%% unfiltered cell figures

figure, hist(I1averageintensity), title('histogram of average fluorescent intensity in each cell');
    xlabel('mean fluorescent intensity (AU)');
    ylabel('number of cells');
    axis([0,(8*10^4),0,500]);
    
    if save == 1
         saveas(gcf, fullfile(fpath,'histogram of average HCR-FISH fluorescent intensity in identified cells (without watershed)'), 'tif');
    end

     [histFreq, histXout] = hist(I1averageintensity, numOfBins);
    figure;
    bar(histXout, histFreq/sum(histFreq)*100), title('percent of cells in average fluorescent intensity bins');
    xlabel('mean fluorescent intensity (AU)');
    ylabel('percent of cells');
    axis([0,(8*10^4),0,100]);

    if save == 1
        savease(gcf, fullfile(fpath, 'percent of cells in average fluorescent intensity bins'), 'tif');
    end
    
    
    figure, hist(I1sumintensity), title('histogram of sum fluorescent intensity in each cell');
    xlabel('sum fluorescent intensity in each cell (AU)');
    ylabel('number of cells');
    axis([0,(5*10^6),0,500]);
    
    if save == 1
         saveas(gcf, fullfile(fpath,'histogram of sum HCR-FISH fluorescent intensity in identified cells (without watershed)'), 'tif');
    end
    
    
    figure, errorbar(I1averageintensity, standardDev, 'rx'), title('average fluorescence and st dev per cell');
    xlabel('cell number');
    ylabel('mean fluorescent intensity (AU) with standard deviation');
    axis([0,400,(-5000),(8*10^4)]);
    
    if save == 1
         saveas(gcf, fullfile(fpath,'average HCR-FISH fluorescent intensity and standard dev'), 'tif');
    end
    
figure, scatter((1:number_cells),I1averageintensity, 'rx'), title('average fluorescence per cell'); %added on 5/12/2016
    xlabel('cell number');
    ylabel('mean fluorescent intensity (AU)');
    axis([0,400,0,(8*10^4)]);
    
    if save == 1
        saveas(gcf, fullfile(fpath, 'average HCR-FISH fluorescent intensity'), 'tif'); 
    end
    
       figure, errorbar(I1sumintensity, standardDev, 'rx'), title('sum fluorescence and st dev');
    xlabel('cell number');
    ylabel('sum fluorescent intensity (AU) with standard deviation of mean');
    axis([0,400,0,(5*10^6)]);
    
    if save == 1
         saveas(gcf, fullfile(fpath,'sum HCR-FISH fluorescent intensity and standard dev'), 'tif');
    end
    
    
%% filtered cells figures

    
     figure, hist(I1averageintensityfiltered), title('histogram of average fluorescent intensity in each cell- filtered');
    xlabel('mean fluorescent intensity (AU)');
    ylabel('number of cells');
    axis([0,(8*10^4),0,400]);
    
    if save == 1
         saveas(gcf, fullfile(fpath,'histogram of average HCR-FISH fluorescent intensity in identified cells (without watershed) filtered'), 'tif');
    end
    
    [histFreq, histXout] = hist(I1averageintensityfiltered, numOfBins);
    figure;
    bar(histXout, histFreq/sum(histFreq)*100), title('percent of cells in average fluorescent intensity bins- filtered');
    xlabel('mean fluorescent intensity (AU)');
    ylabel('percent of cells');
    axis([0,(8*10^4),0,100]);

    if save == 1
        savease(gcf, fullfile(fpath, 'percent of cells in average fluorescent intensity bins- filtered'), 'tif');
    end
    
    figure, hist(I1sumintensityfiltered), title('histogram of sum fluorescent intensity in each cell- filtered');
    xlabel('sum fluorescent intensity in each cell (AU)');
    ylabel('number of cells');
    axis([0,(5*10^6),0,400]);
    
    if save == 1
         saveas(gcf, fullfile(fpath,'histogram of sum HCR-FISH fluorescent intensity in identified cells (without watershed) filtered'), 'tif');
    end
    
        figure, errorbar(I1averageintensityfiltered, standardDevfiltered, 'rx'), title('average fluorescence and st dev per cell- filtered');
    xlabel('cell number');
    ylabel('mean fluorescent intensity (AU) with standard deviation');
    axis([0,400,(-5000),(8*10^4)]);
    
    if save == 1
         saveas(gcf, fullfile(fpath,'average HCR-FISH fluorescent intensity and standard dev filtered'), 'tif');
    end
    
    figure, scatter((1:number_cells_filtered),I1averageintensityfiltered, 'rx'), title('average HCR-FISH fluorescent intensity in each cell'); %added on 5/12/2016
    xlabel('cell number');
    ylabel('mean fluorescent intensity (AU)');
    axis([0,400,0,(8*10^4)]);
    
    if save == 1
        saveas(gcf, fullfile(fpath, 'average HCR-FISH fluorescent intensity filtered'), 'tif'); 
    end
    

    figure, errorbar(I1sumintensityfiltered, standardDevfiltered, 'rx'), title('sum HCR-FISH fluorescent intensity and standard dev filtered');
    xlabel('cell number');
    ylabel('sum fluorescent intensity (AU) with standard deviation of mean');
    axis([0,400,0,(5*10^6)]);
    
    if save == 1
         saveas(gcf, fullfile(fpath,'sum HCR-FISH fluorescent intensity and standard dev filtered'), 'tif');
    end
    
    
%% multiple fields figures

 figure, hist(fieldsaverageintensity), title('histogram of average fluorescent intensity in each cell- multiple filtered fields');
    xlabel('mean fluorescent intensity (AU)');
    ylabel('number of cells');
    axis([0,(5*10^4), 0, inf]); %set x-axis but not y-axis... it will default to max # of cells in that bin
    %axis([0,(5*10^4),0,500]);
    
    if save == 1
         saveas(gcf, fullfile(fpath,'histogram of average HCR-FISH fluorescent intensity in identified cells- multiple fields (without watershed) filtered'), 'tif');
    end

        [histFreq, histXout] = hist(fieldsaverageintensity, numOfBins);
    figure;
    bar(histXout, histFreq/sum(histFreq)*100), title('percent of cells in average fluorescent intensity bins- multiple filtered fields');
    xlabel('mean fluorescent intensity (AU)');
    ylabel('percent of cells');
    axis([0,(8*10^4),0,100]);

    if save == 1
        saveas(gcf, fullfile(fpath, 'percent of cells in average fluorescent intensity bins- multiple fields filtered'), 'tif');
    end
    
        figure, hist(fieldssumintensity), title('histogram of sum fluorescent intensity in each cell- multiple filtered fields');
    xlabel('sum fluorescent intensity in each cell (AU)');
    ylabel('number of cells');
    axis([0,(5*10^6),0,inf]);
    %axis([0,(5*10^6),0,200]);
    
    if save == 1
         saveas(gcf, fullfile(fpath,'histogram of sum HCR-FISH fluorescent intensity in identified cells- multiple fields (without watershed)filtered'), 'tif');
    end
    
     figure, errorbar(fieldsaverageintensity, fieldsstandardDev, 'rx'), title('average fluorescence and st dev per cell- multiple filtered fields');
    xlabel('cell number');
    ylabel('mean fluorescent intensity (AU) with standard deviation');
    axis([0,inf,(-5000),(8*10^4)]); 
    %axis([0,500,(-2000),(6*10^4)]);
    
    if save == 1
         saveas(gcf, fullfile(fpath,'average HCR-FISH fluorescent intensity and standard dev- multiple fields filtered'), 'tif');
    end
    
       number_cells_fields = length(fieldssummary(:,1)); 
    figure, scatter((1:number_cells_fields),fieldsaverageintensity, 'rx'), title('average fluorescent intensity in each cell- multiple fields'); %added on 5/12/2016
    xlabel('cell number');
    ylabel('mean fluorescent intensity (AU)');
    axis([0,inf,0,(8*10^4)]);
    
    if save == 1
        saveas(gcf, fullfile(fpath, 'average HCR-FISH fluorescent intensity multiple filtered fields'), 'tif'); 
    end
    
        figure, errorbar(fieldssumintensity, fieldsstandardDev, 'rx'), title('sum HCR-FISH fluorescent intensity and st dev- filtered fields');
    xlabel('cell number');
    ylabel('sum fluorescent intensity (AU) with standard deviation of mean');
    axis([0,inf,0,(5*10^6)]);
    %axis([0,500,0,(1*10^7)]);
    
    if save == 1
         saveas(gcf, fullfile(fpath,'sum HCR-FISH fluorescent intensity and standard dev- multiple fields filtered'), 'tif');
    end 
    
%% multiple field summary

figure, bar(meanareas), title('mean cell area for each field filtered');
    xlabel('cell field');
    ylabel('mean cell area (pixels)');
    
    if save == 1
         saveas(gcf, fullfile(fpath,'mean cell area for each field filtered'), 'tif');
    end
    
    figure, bar(averagemeanintensities), title('average mean intensities of cells in each field filtered ');
    xlabel('cell field');
    ylabel('average of mean cell intensities (AU)');
    
    if save == 1
         saveas(gcf, fullfile(fpath,'average mean intensities of cells in each field filtered'), 'tif');
    end
    
    figure, bar(averagesumintensities), title('average sum intensities of cells in each field filtered');
    xlabel('cell field');
    ylabel('average of sum cell intensity (AU)');
    
    if save == 1
         saveas(gcf, fullfile(fpath,'average sum intensities of cells in each field filtered'), 'tif');
    end
    
    figure, bar(averagestandarddevs), title('average standard dev of fluor intensity of cells in each field filtered');
    xlabel('cell field');
    ylabel('average of standard dev of cell intensity');
    
    if save == 1
         saveas(gcf, fullfile(fpath,'average standard dev of fluor intensity of cells in each field filtered'), 'tif');
    end

    