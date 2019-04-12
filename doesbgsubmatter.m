%% initiation code

clc
clear

save = 1; %set to 1 if you want to save figures
if save == 1; %need to set foldername to save to 
mkdir('Sintensity_test1') %replace foldername with desired output folder
end
fpath = '\\Client\C$\1MATLABanalysis\EfacDAPILacZ100x\Sintensity_test1'; %edit foldername here too!

%% load images

[Iraw, fullpathraw] = loadAndShowImage(1,1);
[I, fullpath] = loadAndShowImage(1,1);
[I1raw, fullpath1raw] = loadAndShowImage(1,1);
[I1, fullpath1] = loadAndShowImage(1,1);

%% plot profile 

[x,y] = ginput(2); %get clicks for two points from which a line will be drawn
xline = [x(1), x(2)];
yline = [y(1), y(1)]; %ignore noise along y axis

figure
improfile(I, xline, yline), grid on, title ('arbitrary intensity profile, background subtracted image');
saveas(gcf, fullfile(fpath,'1 intensity Hoechst'), 'tif');

figure
improfile(Iraw, xline, yline), grid on, title ('arbitrary intensity profile, raw image');
saveas(gcf, fullfile(fpath,'2 intensity Hoechst raw'), 'tif');

figure
improfile(I, xline, yline), grid on, title ('arbitrary intensity profile, background subtracted image 1');
saveas(gcf, fullfile(fpath,'3 intensity lacZ'), 'tif');

figure
improfile(Iraw, xline, yline), grid on, title ('arbitrary intensity profile, raw image 1');
saveas(gcf, fullfile(fpath,'4 intensity lacZ raw'), 'tif');
