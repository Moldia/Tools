%% Tile one image series
%  Xiaoyan, 2017


%% parameters
imagename_prefix = 'Anne_B_032818_c';
imagename_suffix = '_ORG.tif';
series_num = 2;
tilesize = 4000;

%% read and pad the image
for s = 1:series_num
    I = imread([imagename_prefix, num2str(s), imagename_suffix]);
    outputdir = ['Tiled_',...
        strtok([imagename_prefix, num2str(s), imagename_suffix],'.')];
    tileimage(I, tilesize, outputdir);
end