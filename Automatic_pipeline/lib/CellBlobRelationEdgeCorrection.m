% correct cell objects cut by tiling line
% re-relate blob and parent
% Xiaoyan, 2017


%% input
tilepos = getcsvtilepos('I:\MIPs_ms120\CellProfiler\Tiled_wref.csv');
blobfile = 'I:\MIPs_ms120\CellProfiler\Decoding\QT_0.45_0.0001_details.csv';

%% stitching and edge correction
tilesize = 2000;
ntilesX = max(tilepos(:,2)/tilesize) + 1;
ntilesY = max(tilepos(:,3)/tilesize) + 1;

ndigits = 0;

mkdir('Stitched');
mkdir('ParentCell')
mkdir('CellBlobs');

% relabel tile object images outputed from CP
Ilabeled = cell(ntilesY,ntilesX);
Ioutlines = cell(ntilesY,ntilesX);

for i = 1:ntilesY
    fprintf('%.2f percent read.\n', (i-1)/ntilesY*100);
    parfor j = 1:ntilesX
        t = ntilesX*(i-1)+j;
        I = imread(['Segmentation\tile', num2str(t), '_Nuclei.tif']);
        Ivis = imread(['Outlines\tile', num2str(t), '.png']);
        I = relabelimg(I);
        Ilabeled(i,j) = {I};
        Ioutlines(i,j) = {Ivis};
    end
    
end  
clear I Ivis
fprintf('100 percent finished.\n');

% stitched outline images for visualization
resizef = .5;
Ivis = resizestitch([ntilesX,ntilesY], tilesize, resizef, Ioutlines);

%% relate blobs to parent cell (if not done in CellProfiler)
Icell = cell(ntilesY,ntilesX);
for i = 1:ntilesY
    fprintf('%.2f percent read.\n', (i-1)/ntilesY*100);
    parfor j = 1:ntilesX
        t = ntilesX*(i-1)+j;
        I = imread(['Segmentation\tile', num2str(t), '_Cell.tif']);
        Icell(i,j) = {relabelimg(I)};
    end  
end  
clear I
fprintf('100 percent finished.\n');

[~, pos] = getinsitudata(blobfile);
pos = correctcoord(pos, 1);     % difference between zero indexing (Python) and non-zero indexing
cellid = zeros(length(pos),1);
for i = 1:ntilesY
    for j = 1:ntilesX
        in = readsinsqr(pos, tilesize*[j-1 j i-1 i]+.5);
        postmp = bsxfun(@minus, pos(in,:), tilesize*[j-1 i-1]);
        postmp = readsinroi(postmp, Icell{i,j});
        cellid(in) = postmp;
    end
end
writeblobwcell(blobfile, cellid, 'QT0.35');

%% cell/nuclei properties
cellprop = importdata('Nuclei.csv');     % needs area, centroid position
% cellprop: image number, object number, tile number, area, x, y
cellprop = [cellfun(@str2num, cellprop.textdata(2:end,1)),...
    cellfun(@str2num, cellprop.textdata(2:end,2)),...
    cellprop.data(:,[3,4:6])];

% cell global position
if size(cellprop,2) > 3
    for i = 1:max(cellprop(:,3))
        try
            cellprop(cellprop(:,3)==i, end-1:end) = ...
                bsxfun(@plus, cellprop(cellprop(:,3)==i,end-1:end), tilepos(i,2:3));
        end
    end
else
    warning('Too few columns detected in cell property file.') 
end

% cellprop(:,[1,2]) = cellprop(:,[2,1]);

%% find objects cut by tiling lines
[CellCorrectionTable, CellProps, Ioutlines] = correct_edge_objects...
    (Ilabeled, Ioutlines, cellprop, tilepos, Ivis, resizef);
save('Stitched\CellLookupTable.mat', 'CellCorrectionTable', 'CellProps');

% remake the resized outline image
Ivis = resizestitch([ntilesX,ntilesY], tilesize, resizef, Ioutlines);
imwrite(Ivis, ['Stitched\Outlines_' num2str(resizef*100) '%.jpg']);
IvisHD = resizestitch([ntilesX,ntilesY], tilesize, 1, Ioutlines);
imwrite(IvisHD, 'Stitched\Outlines_100%.jpg');

%% stitch full section DAPI label image and expand
% IlabelCorrected = Ilabeled;
% parfor i = 1:length(IlabelCorrected(:))
%     if ~isempty(CellCorrectionTable{i})
%         for j = fliplr(CellCorrectionTable{i}(:,1)')
%             IlabelCorrected{i}(IlabelCorrected{i} == j) = ...
%                 CellCorrectionTable{i}(CellCorrectionTable{i}(:,1)==j,2);
%         end
%     end
% end
% LabelDapiWhole = resizestitch([ntilesX,ntilesY], tilesize, resizef, IlabelCorrected);
% 
% [D, idx] = bwdist(LabelDapiWhole);
% D = uint32(D<=20);
% L = reshape(LabelDapiWhole(idx(:)), size(LabelDapiWhole));
% LabelCellWhole = D.*L;
% 
% IOutlines = LabelCellWhole ~= imerode(LabelCellWhole, strel('disk', 2));
% 
% imwrite(IOutlines, 'Stitched\Outlines.jpg');
% save('Stitched\LabelImages', 'LabelCellWhole', 'LabelDapiWhole');

% clear Ilabeled Ioutlines IvisHD

%% correct cell properties (area, centroid position)
%  properties should be from expanded objects ("real parents")
figure,imshow(Ivis)
hold on

% renumber cells
cellrenumber = renumberedgeobj...
    ([ntilesX,ntilesY], CellCorrectionTable, CellProps, 2);
[uniCell, ~, idxCell] = unique(cellrenumber, 'stable');
countCell = hist(idxCell, 1:length(uniCell));
fid = fopen('Stitched\UniqueCells.txt', 'w');
fprintf(fid, '%d\n', uniCell);
fclose(fid);

CellPropsRenum = [uniCell(countCell==1),...
    CellProps(ismember(idxCell, find(countCell==1)),3:6)];
plot(CellPropsRenum(:,end-1)*resizef, CellPropsRenum(:,end)*resizef, '.');

% merge cells and update properties
edgecells = uniCell(countCell~=1);
for i = 1:length(edgecells)
    if edgecells(i) > 0
        edgecellProps = CellProps(cellrenumber==edgecells(i),3:6);
        % reassign tile number according to biggest part and sum up area
        [~, tilere] = max(edgecellProps(:,2));
        tile = edgecellProps(:,1);
        tilere = tile(tilere);
        atotal = sum(edgecellProps(:,2),1);
        plot(edgecellProps(:,3)*resizef,edgecellProps(:,4)*resizef,'-');
        % recalculate center of mass
        pos = edgecellProps(:,3:4).*repmat(edgecellProps(:,2)/atotal, 1, 2);
        try
            pos = sum(pos,1);
        catch
            pos = [nan, nan];
        end
        plot(pos(:,1)*resizef,pos(:,2)*resizef,'o');
        CellPropsRenum = [CellPropsRenum; [edgecells(i), tilere, atotal, pos]];
    end
end

[~,idxsort] = sort(CellPropsRenum(:,1));
CellPropsRenum = CellPropsRenum(idxsort,:);
fid = fopen('Stitched\Cells.csv', 'w');
fprintf(fid, 'CellID,metadata_position,area,global_x_pos,global_y_pos\n');
fprintf(fid, '%d,%d,%d,%d,%d\n', reshape(CellPropsRenum',[],1));
fclose(fid);
save('Stitched\CellLookupTable.mat', 'CellPropsRenum', '-append');

%% parent cell based on expanded full section DAPI label
% [~, pos] = getinsitudata(blobfile);
% pos = correctcoord(pos, 1);     % difference between zero indexing (Python) and non-zero indexing
% parentcellNew = readsinroi(pos, LabelCellWhole);
% 
% writeblobwcell(blobfile, parentcellNew, 'beforeQT');

%% blob cell relation     
parentcell = csvread('ParentCell\QT0.35_details_wCell.csv', 1, 3);
columnCell = 7;

% % original tile size different cell segmentation tile size
% tileid = ceil(pos/2000);
% tileid = (tileid(:,2)-1)*ntilesX + tileid(:,1);
% parentcell(:,3) = tileid;

parentcellNew = renumberedgeobj...
    ([ntilesX,ntilesY], CellCorrectionTable, parentcell, columnCell);

towrite = [(1:size(parentcell,1))', parentcell(:,3), parentcellNew]';
fid = fopen('ParentCell\ParentCell_QT0.35Blobs.csv', 'w');
fprintf(fid,'blob_num, metadata_position, Parent_Cell\n');
fprintf(fid, '%d,%d,%d\n', towrite(:));
fclose(fid);                                                                                                                                    

%% label child blobs and get single cell profile
% before QT
% parent = findparentcell('..\CP_170122_Seq\Decoding\beforeQT.mat',...
%     'blob_allbt', parentcellNew);
parent = parentcellNew;
csvwrite('ParentCell\ParentCell_QT0.35.csv', parent);
writeblobwcell(blobfile, parent, 'QT0.35_b124');
childblobs('ParentCell\QT0.35_details_wCell.csv',...
    parent, CellPropsRenum, 'QT0.35_b124');

