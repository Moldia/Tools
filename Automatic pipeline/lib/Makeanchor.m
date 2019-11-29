PATH='J:\MIPs_ms120\190826_NewStitch\Export\AlignedImages\'
FILE_PREFIX='Cyc1m1_'
FILE_SUFIX='_ORG.tif'




for cycle = 1:5
for tile = 1:3
cycle=1
tile=1
FILE_SUFIX=['_ORG_align.tif']    
FILE_PREFIX=['Cyc',num2str(cycle),'m',num2str(tile),'_']
create_anchor (PATH,FILE_PREFIX, FILE_SUFIX)

end 
end 


function create_anchor (PATH,FILE_PREFIX, FILE_SUFIX)

flo1 = imread([PATH, FILE_PREFIX,'2',FILE_SUFIX]);
flo2 = imread([PATH, FILE_PREFIX,'3',FILE_SUFIX]);
flo3 = imread([PATH, FILE_PREFIX,'4',FILE_SUFIX]);
flo4 = imread([PATH, FILE_PREFIX,'5',FILE_SUFIX]);
 
 %% Tophat
 flo_top1= imtophat(flo1,strel('disk', 2));
 flo_top2= imtophat(flo2,strel('disk', 2));
 flo_top3= imtophat(flo3,strel('disk', 2));
 flo_top4= imtophat(flo4,strel('disk', 2));
 th=100;
 
v=0
 BW=flo_top1>th;
 s=regionprops(BW);
 coordinates2 = [s(:).Centroid];
 sizec=size(coordinates2);
 ematrix = reshape(coordinates2,2,sizec(2)/2);
 ematrix=transpose(ematrix);
detection1=zeros(size(flo_top1));
size1=size(ematrix);
for i=1:size1(1)
   detection1(round(ematrix(i,1)),round(ematrix(i,2)))=100; 
end
 

 
v=1


%TOP2
  BW=flo_top2>th;
 s=regionprops(BW);
 coordinates2 = [s(:).Centroid];
 sizec=size(coordinates2);
 ematrix = reshape(coordinates2,2,sizec(2)/2);
 ematrix=transpose(ematrix);
detection2=zeros(size(flo_top2));
size2=size(ematrix);
for i=1:size2(1)
   detection2(round(ematrix(i,1)),round(ematrix(i,2)))=100; 
end
 
 
v=2
%top3
  BW=flo_top3>th;
 s=regionprops(BW);
 coordinates3 = [s(:).Centroid];
 sizec=size(coordinates3);
 ematrix = reshape(coordinates3,2,sizec(2)/2);
 ematrix=transpose(ematrix);
detection3=zeros(size(flo_top3));
size3=size(ematrix);
for i=1:size3(1)
   detection3(round(ematrix(i,1)),round(ematrix(i,2)))=100; 
end
 
v=3

%top3
   BW=flo_top4>th;
 s=regionprops(BW);
 coordinates4 = [s(:).Centroid];
 sizec=size(coordinates4);
 ematrix = reshape(coordinates4,2,sizec(2)/2);
 ematrix=transpose(ematrix);
detection4=zeros(size(flo_top4));
size4=size(ematrix);
for i=1:size4(1)
   detection4(round(ematrix(i,1)),round(ematrix(i,2)))=100; 
end

v=4

anchor= detection1+detection2+detection3+detection4;
anchor_filt= imtophat(anchor,strel('disk', 2));
%contrastAdjusted = imadjust(anchor_filt);

imwrite(anchor_filt, [PATH, FILE_PREFIX,'6',FILE_SUFIX,'_nw2']);
end
 
 