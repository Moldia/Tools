% kernel density estimation of all genes
% Xiaoyan, 2017

clear;
close all;

%% parameters
decoded_file = 'D:\Salamander project\Ch030\CP_170418\Decoding\QT_0.55_0.0001_details_ROI.csv';
image = 'D:\Salamander project\Ch030\Preprocessing\Stitched\Ch030_170310_Pw_5dpa_SBL2_CX_c1_stitched.tif';   % important for size
scale = 1;      % image scale
bandwidth = 30; % in original scale


%%
% transcripts
[name, pos] = getinsitudata(decoded_file);
[name, pos] = removereads(name, 'NNNN', pos);

% unique transcripts
[uNames, ~, idxName] = unique(name);
cNames = hist(idxName, 1:length(uNames));

% % remove reads with low counts
% [name, pos] = removereads(name, uNames(cNames<200), pos);
% [uNames, ~, idxName] = unique(name);
% cNames = hist(idxName, 1:length(uNames));

% density estimation plot
f = 0;
figure;
set(gcf, 'units', 'normalized', 'position', [.05+.02*f .1-.02*f .8 .8]);
for i = 1:length(uNames)
    if mod(i-1,12)==0 && i~=1
        drawnow;
        f = f+1;
        saveas(gcf, ['Density_', num2str(bandwidth), '_', num2str(f), '.png']);
        disp([num2str(12*f) ' transcripts are finished.']);
        
        clf;
        set(gcf, 'units', 'normalized', 'position', [.05+.02*f .1-.02*f .8 .8]);
    end
    density = gene_kde(name, pos, uNames{i}, bandwidth, image, scale);
    
    subplot(3, 4, i-f*12);
    imshow(density, []);
    colormap(gca, parula);
    title([uNames{i}, ' (', num2str(cNames(i)), ')']);
end
disp('All transcripts are finished.')
saveas(gcf, ['Density_', num2str(bandwidth), '_', num2str(f+1), '.png']);

