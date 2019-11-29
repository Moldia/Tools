% kernel density estimation of all genes
% Xiaoyan, 2017

clear;
close all;

%% parameters
decoded_file = 'K:\Organoid_method_development\Organoids_erik\Organoid3_Output\1\QT_0.6_details_noNNNN.csv';
image = 'K:\Organoid_method_development\Organoids_erik\Organoid3_Output\1\Sample 1_Base 3_resize20_c1.jpg';   % important for size
scale = 0.2;      % image scale
bandwidth = 100;   % in original scale


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

