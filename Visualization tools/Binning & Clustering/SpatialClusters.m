% DBSCAN (Density-based spatial clustering of applications with noise)
% Xiaoyan, 2017


%% modify here
[name, pos] = getinsitudata_f('E:\PROOOJECTS\11_PDGFR\K\CP_150708\Decoding\QT_0.4_0.001_details.csv');
[uNames, ~, iName] = unique(name);
background = imread('E:\PROOOJECTS\11_PDGFR\K\IHC\IHC_mergedRGB_20%.tif');
scale = 0.2;

name_density = 'CDH1';
max_dist = 700;

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

