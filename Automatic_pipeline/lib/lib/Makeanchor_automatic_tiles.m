PATH='J:\MIPs_ms120\sergio\'
OUTPATH='J:\MIPs_ms120\sergio\withanchor\'
FILE_PREFIX='190917_msSBH_02_rd6cy1(1)_t1.tif'
FILE_SUFIX='_ORG.tif'
filenames={
    'HCA_09_ms_120_c1-MIP_sergio\HCA_09_ms_120_c1-MIP_sergio_c',
    'HCA_09_ms_120_c2-MIP\HCA_09_ms_120_c2-MIP_c',
    'HCA_09_ms_120_c3-MIP\HCA_09_ms_120_c3-MIP_c',
    'HCA_09_ms_120_c4-MIP\HCA_09_ms_120_c4-MIP_c',
    'HCA_09_ms_120_c5-MIP\HCA_09_ms_120_c5-MIP_c',
}

for tile = 1:9
for cyc= 1:size(filenames)   

FILE_SUFIX=['m',num2str(tile),'_ORG.tif']    
FILE_PREFIX=strcat(filenames(cyc))
class(FILE_PREFIX)
FILE_PREFIX=char(FILE_PREFIX)
create_anchor (PATH,FILE_PREFIX, FILE_SUFIX,tile,cyc,OUTPATH)

end 
end


function create_anchor (PATH,FILE_PREFIX, FILE_SUFIX,tile,cyc,OUTPATH)
flodap = imread(strcat(PATH, FILE_PREFIX,'1',FILE_SUFIX));
flo1 = imread(strcat(PATH, FILE_PREFIX,'2',FILE_SUFIX));
flo2 = imread(strcat(PATH, FILE_PREFIX,'3',FILE_SUFIX));
flo3 = imread(strcat(PATH, FILE_PREFIX,'4',FILE_SUFIX));
flo4 = imread(strcat(PATH, FILE_PREFIX,'5',FILE_SUFIX));
 
 %% Tophat
 flo_top1= imtophat(flo1,strel('disk', 2));
 flo_top2= imtophat(flo2,strel('disk', 2));
 flo_top3= imtophat(flo3,strel('disk', 2));
 flo_top4= imtophat(flo4,strel('disk', 2));
 th=100;
 
% v=0
%  BW=flo_top1>th;
%  s=regionprops(BW);
%  coordinates2 = [s(:).Centroid];
%  sizec=size(coordinates2);
%  ematrix = reshape(coordinates2,2,sizec(2)/2);
%  ematrix=transpose(ematrix);
% detection1=zeros(size(flo_top1));
% size1=size(ematrix);
% for i=1:size1(1)
%    detection1(round(ematrix(i,1)),round(ematrix(i,2)))=100; 
% end
%  
% 
%  
% v=1
% 
% 
% %TOP2
%   BW=flo_top2>th;
%  s=regionprops(BW);
%  coordinates2 = [s(:).Centroid];
%  sizec=size(coordinates2);
%  ematrix = reshape(coordinates2,2,sizec(2)/2);
%  ematrix=transpose(ematrix);
% detection2=zeros(size(flo_top2));
% size2=size(ematrix);
% for i=1:size2(1)
%    detection2(round(ematrix(i,1)),round(ematrix(i,2)))=100; 
% end
%  
%  
% v=2
% %top3
%   BW=flo_top3>th;
%  s=regionprops(BW);
%  coordinates3 = [s(:).Centroid];
%  sizec=size(coordinates3);
%  ematrix = reshape(coordinates3,2,sizec(2)/2);
%  ematrix=transpose(ematrix);
% detection3=zeros(size(flo_top3));
% size3=size(ematrix);
% for i=1:size3(1)
%    detection3(round(ematrix(i,1)),round(ematrix(i,2)))=100; 
% end
%  
% v=3
% 
% %top3
%    BW=flo_top4>th;
%  s=regionprops(BW);
%  coordinates4 = [s(:).Centroid];
%  sizec=size(coordinates4);
%  ematrix = reshape(coordinates4,2,sizec(2)/2);
%  ematrix=transpose(ematrix);
% detection4=zeros(size(flo_top4));
% size4=size(ematrix);
% for i=1:size4(1)
%    detection4(round(ematrix(i,1)),round(ematrix(i,2)))=100; 
% end
% 
% v=4
% 
% 
% anchor= detection1+detection2+detection3+detection4;

anchor= flo_top1+flo_top2+flo_top3+flo_top4;
anchor_filt= imtophat(anchor,strel('disk', 2));
%contrastAdjusted = imadjust(anchor_filt);
imwrite(flodap,strcat(OUTPATH,'Tile',num2str(tile),'_cyc_',num2str(cyc),'_CH_1',FILE_SUFIX))
imwrite(flo1,strcat(OUTPATH,'Tile',num2str(tile),'_cyc_',num2str(cyc),'_CH_2',FILE_SUFIX))
imwrite(flo2,strcat(OUTPATH,'Tile',num2str(tile),'_cyc_',num2str(cyc),'_CH_3',FILE_SUFIX))
imwrite(flo3,strcat(OUTPATH,'Tile',num2str(tile),'_cyc_',num2str(cyc),'_CH_4',FILE_SUFIX))
imwrite(flo4,strcat(OUTPATH,'Tile',num2str(tile),'_cyc_',num2str(cyc),'_CH_5',FILE_SUFIX))
imwrite(anchor_filt, strcat(OUTPATH,'Tile',num2str(tile),'_cyc_',num2str(cyc),'_CH_6',FILE_SUFIX));
end
 
 