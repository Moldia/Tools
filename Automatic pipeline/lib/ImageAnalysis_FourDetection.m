%% Image analysis pipeline for four detection
%  Channels: c1=Cy3, c2=TRed, c3=Cy5, c4=AF750
%  Tailored for lung cancer sections (Ephrin project)
%  Xiaoyan, 2015-1-4


%%
clear
% cd('E:\PROOOJECTS\9_Lung_Multiplex\11692\11692_2_s1_m2');
format compact

%% parameters
Thtop = [.005 .004 .002 .007];
Thback = [.009 .0065 .004 .015];
Rback = [1.5 4 1 2];
    
%% loading and padding
Inuclei = imread('ZENout\11692_2-MIP_s2c1_ORG.tif');
Isize = size(Inuclei);
tileinum = ceil(Isize(1)/2000);
tilejnum = ceil(Isize(2)/2000);
padi = tileinum*2000 - Isize(1);
padj = tilejnum*2000 - Isize(2);
% pad
Inuclei = [Inuclei,zeros(Isize(1),padj);zeros(padi,Isize(2)+padj)];
Itemp = zeros(size(Inuclei));    % preallocate 

% load images
I = [];
for i = 1:4
    Itemp = imread(['ZENout\11692_2-MIP_s2c' num2str(i+1) '_ORG.tif']);
    % channel alignment
    if i == 1
        Itemp = [zeros(1,size(Itemp,2));Itemp(1:end-1,2:end),zeros(size(Itemp,1)-1,1)];
    end
    Itemp = [Itemp,zeros(Isize(1),padj);zeros(padi,Isize(2)+padj)]; % pad
    Itemp = double(Itemp)/65525;
    field = ['c' num2str(i)];
    I.(field) = Itemp;
end
disp('image loading finished');

%% image analysis
tic
clf
Count = zeros(tileinum*tilejnum,6);
Pos = cell(1,4);
counter = 0;
mkdir('Segmentation');
for i = 1:tileinum
    for j = 1:tilejnum
        % tissue region
        Inuc = Inuclei(2000*(i-1)+1:2000*i,2000*(j-1)+1:2000*j);
        Imask = im2bw(Inuc,graythresh(Inuclei));
        Imask = imdilate(Imask,strel('disk',50));
        % imshow(Imask);
        
        Itile = cat(3,...
            I.c1(2000*(i-1)+1:2000*i,2000*(j-1)+1:2000*j),...
            I.c2(2000*(i-1)+1:2000*i,2000*(j-1)+1:2000*j),...
            I.c3(2000*(i-1)+1:2000*i,2000*(j-1)+1:2000*j),...
            I.c4(2000*(i-1)+1:2000*i,2000*(j-1)+1:2000*j));
        
        % segmentation
        I_thresh = [];
        I_threshdil = [];
        I_comp = [];
        for c = 1:4
            Itemp = Itile(:,:,c);
            % Itemp = Itemp.*Imask;
            
            % filters
            Iback = imclose(Itemp,strel('disk',10));
            Iback = imdilate(Iback,strel('disk',4));
            Itop = imtophat(Itemp,strel('disk',5));
            
            % thresholding and size filtering
            Ithreshraw = im2bw(Itop,Thtop(c));
            Ithresh = imfill(Ithreshraw,'holes');
            Ithresh = bwareaopen(Ithresh,2);
            big = bwareaopen(Ithresh,300);
            Ithresh = Ithresh - big;
            Ithresh(Ithresh<0) = 0;
            
            % signals under autofluorescence
            Icomp = zeros(size(Itemp));
            switch c
                case 1
                    Icomp(Itemp>=Itile(:,:,2)*Rback(c)) = 1;
                case 2
                    Icomp(Itemp>=Itile(:,:,1)*Rback(c)) = 1;
                case 3
                    Icomp(Itemp>=Itile(:,:,1)*Rback(c) & Itemp>=Itile(:,:,2) & Itemp>=Itile(:,:,4)*0.5) = 1;
                case 4
                    Icomp(Itemp>=Itile(:,:,1)*Rback(c) & Itemp>=Itile(:,:,2) & Itemp>=Itile(:,:,3)*4) = 1;
            end
            Icomp = bwareaopen(Icomp,5);
            
            % strong background structures
            Ithreshback = im2bw(Iback,Thback(c));
            Ithreshback = Ithreshback - Icomp;
            Ithreshback(Ithreshback<0) = 0;
            Ithreshback = bwareaopen(Ithreshback,150);
            Ithresh = Ithresh - Ithreshback;
            Ithresh(Ithresh<0) = 0;
            
            % remove signals outside tissue region (alternative)
            Ithresh = Ithresh.*Imask;
            
            Ithreshd = imdilate(Ithreshraw,strel('disk',5));
            
            % imshow(cat(3,Itemp*50,Ithresh,Ithresh))
            % imshow(cat(3,Itemp*50,Ithreshback,Ithreshback))
            % imshow(cat(3,Itemp*50,Icomp,Icomp))
            
            I_thresh = cat(3,I_thresh,Ithresh);
            I_threshdil = cat(3,I_threshdil,Ithreshd);
            I_comp = cat(3,I_comp,Icomp);
        end
        
        % signals appearing in multiple channels
        Imulti = mean(I_threshdil,3);
        Imulti(Imulti>.25) = 1;
        Imulti(Imulti==.25) = 0;    
        % imshow(cat(3,Imulti,zeros(size(Itemp)),zeros(size(Itemp))) + cat(3,double(Inuc)/65525,double(Inuc)/65525,double(Inuc)/65525))
        
        % remove signals appearing in multiple channels
        count = [];
        for c = 1:4
            % signals found in multiple channels (dilated segmetation)
            Itemp = imreconstruct(Imulti,I_threshdil(:,:,c));
            % signals found in multiple channels (original segmentation)
            Itemp = imreconstruct(Itemp,I_thresh(:,:,c));
            % remove the ones should be compensated
            Itemp = Itemp - I_comp(:,:,c);
            Itemp(Itemp<0) = 0;
            % imshow(I_comp(:,:,c))
            
            Itemp = I_thresh(:,:,c) - Itemp;
            Itemp(Itemp<0) = 0;
            Itemp = imfill(Itemp,'holes');
            
            I_thresh(:,:,c) = Itemp;
            
            % count objects
            obj = bwconncomp(Itemp);
            count = [count,obj.NumObjects];
            
            % the first index of an object (approximate position)
            p = cellfun(@(p) p(1),obj.PixelIdxList);
            p = [ceil(p'/2000),mod(p',2000)];
            p(p(:,1)==0) = 2000;
            p(p(:,2)==0) = 1;
            p = p + repmat([j-1,i-1]*2000,size(p,1),1);
            
            % size of each object
            area = cellfun(@length,obj.PixelIdxList);

            Pos{1,c} = [Pos{1,c};p,area'];
        end
        
        % H = figure;
        imshow(cat(3,I_thresh(:,:,2),I_thresh(:,:,1),zeros(size(Itemp))) + ...
            cat(3,I_thresh(:,:,4),zeros(size(Itemp)),I_thresh(:,:,4)) + ...
            cat(3,I_thresh(:,:,3),I_thresh(:,:,3)*.3,zeros(size(Itemp))) + ...
            cat(3,double(Inuc)/65525,double(Inuc)/65525,double(Inuc)/65525))

        counter = counter+1;
        Count(counter,:) = [i,j,count];
        % waitfor(H)
        saveas(gcf,['Segmentation\' num2str(i) '_' num2str(j) '.png'],'png');
        disp([num2str(counter/tileinum/tilejnum*100) '% of images are processed.']);
    end
end
toc
save('ImageAnalysis.mat','Pos','Count','Isize')

%% merge signals of the same color that are very close to each other
% can be a result of re-tiling (one blob separated into two tiles)
for c = 1:4
    NS_idx = rangesearch(Pos{1,c}(:,1:2),Pos{1,c}(:,1:2),30);
    NS_num = cell2mat(cellfun(@size,NS_idx,'UniformOutput',false));
    if ~isempty(NS_num)
        NS_true = find(NS_num(:,2)~=1);
    else
        NS_true = [];
    end
    % modify the pixel size to 0
    for i = 1:length(NS_true)
        if Pos{1,c}(NS_true(i),3) ~= 0
            temp = NS_idx{NS_true(i),1}(2:end);
            Pos{1,c}(temp,3) = 0;
        end
    end
end
Pos_re = Pos;
save('ImageAnalysis.mat','Pos_re','-append')

%% plotting
col = [0 1 0; 1 0 0; 1 .3 0; 1 0 1];
Title = {'EGFR mut' 'A2 wt' 'EGFR wt' 'A2 mut'};
H = figure;
imshow(Inuclei,[0 30000]);
hold on;
Hde = figure;
Hgau = figure;

pos_re = [];
for c = 1:4
    pos = Pos_re{1,c};
    pos = pos(pos(:,3)>=4,:);
    
    pos_re = [pos_re;pos,repmat(c,size(pos,1),1)];
    
    if size(pos,1)>2

        % 2d kernel
        [bandwidth,density,X,Y]=kde2d_X(pos(:,1:2),2^8,[0 0],[Isize(2) Isize(1)],[300 300]);
        
        % gaussian smoothing
        temp = floor(pos(:,1:2)/5);
        temp(temp==0) = 1;
        Itemp = accumarray(fliplr(temp),1,floor(Isize/5));
        fh = fspecial('gaussian',500,50);
        Itemp = imfilter(Itemp,fh);
        
    else
        density = zeros(2^8);
        Itemp = imresize(density,[Isize(1)/5 Isize(2)/5]);
    end
   
    density = imresize(density,[Isize(1)/5 Isize(2)/5]);
    figure(Hde); subplot(2,2,c);
    imshow(density,[0 5E-8]);
    colormap(parula);
    title(Title{c});
    
    figure(Hgau); subplot(2,2,c);
    imshow(Itemp,[1E-5 1E-3]);
    colormap(hot);
    title(Title{c});
    
    % plotting
    figure(H);
    plot(pos(:,1),pos(:,2),'.','color',col(c,:));
end
drawnow;
saveas(Hde,'density_estimation.png','png');
savefig(Hde,'density_estimation.fig');
saveas(H,'Plotted.png','png');
savefig(H,'Plotted.fig');
saveas(Hgau,'gaussian_smoothing.png','png');
savefig(Hgau,'gaussian_smoothing.fig');
save('ImageAnalysis.mat','pos_re','-append')

%% output files
write = [];
for c = 1:4
    write = [write;Pos{1,c},repmat(c,size(Pos{1,c},1),1)];
end
write = write';

fid = fopen('Decoded_original_details.csv','w');
fprintf(fid,'X_pos,Y_pos,px_area,code\n');
fprintf(fid,'%d,%d,%d,%d\n',write(:));
fclose(fid);

fid = fopen('Decoded_details.csv','w');
fprintf(fid,'X_pos,Y_pos,px_area,code\n');
fprintf(fid,'%d,%d,%d,%d\n',pos_re');
fclose(fid);

fid = fopen('Decoded_original_count.csv','w');
fprintf(fid,'tile_i,tile_j,count_1,count_2,count_3,count_4\n');
fprintf(fid,'%d,%d,%d,%d,%d,%d\n',Count');
Total = sum(Count,1);
fprintf(fid,'\n');
fprintf(fid,'sum,,%d,%d,%d,%d',Total(3:end));
fclose(fid);

totalcount = hist(pos_re(:,4),1:4);
fid = fopen('Decoded_count.csv','w');
fprintf(fid,'count_1,count_2,count_3,count_4\n');
fprintf(fid,'%d,%d,%d,%d\n',totalcount);
fclose(fid);