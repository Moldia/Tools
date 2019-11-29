% DBSCAN (Density-based spatial clustering of applications with noise)
% Xiaoyan, 2017


%% modify here
[name, pos] = getinsitudata_f('E:\Whole_organoid_pseudoAnchor_v5\Decoding\QT_0.7_details_noNNNN.csv');
[uNames, ~, iName] = unique(name);
background = imread('E:\Whole_organoid_pseudoAnchor_v5\Align\Base_1_aligned-1.tif');
scale = 1;

name_density = 'Homo sapiens tubulin beta 3 class III (TUBB3)';
max_dist = 100;

%% do not modify

idx_density = find(strcmp(uNames,name_density));

if isempty(idx_density)
    error('No specified transcript detected in the input file');
end
pos_density = pos(iName==idx_density,1:2);
[lab,labc] = dbscan(pos_density,max_dist,20);

pos_density = correctcoordinates_f(pos_density,scale);

figure,imshow(background),hold on;
for i = 1:max(lab)
    plot(pos_density(lab==i,1),pos_density(lab==i,2),'+');
end

