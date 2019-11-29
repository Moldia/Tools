% visualize expression level in tSNE map and in tissue
% press any key to continue after each time the plot is updated
% Xiaoyan, 2018

clear;
close all;

%% modify here
tsne_3d_file = 'J:\NextPipeline\tSNE_3D.csv';
hexbin_file = 'J:\NextPipeline\HexbinClustering_GeneCount.csv';
hexbin_size = 300;
genes_to_show = {''};     % if empty, all genes will be shown one by one

%% do not modify

% load bin count
data = importdata(hexbin_file);

genes = data.colheaders(:,2:end-2);
pos = data.data(:,end-1:end);
counts = data.data(:,2:end-2);

% load tSNE results
tsne_3d = csvread(tsne_3d_file);
tsne_rgb = rgbscale(tsne_3d);

if isempty(genes_to_show)
    genes_to_show = genes;
end

figure; 
for i = 1:numel(genes_to_show)
    Ax = [];
    idx = find(strcmp(genes_to_show{i}, genes));
    if ~isempty(idx)
        col = summer(max(counts(:,idx))+1);

        clf;

        % tSNE colored by expression
        ax = subplot(221); hold on;
        for j = 1:size(pos,1)
            plot3(tsne_3d(j,1), tsne_3d(j,2), tsne_3d(j,3), '.',...
                'color', col(counts(j,idx)+1,:));
        end
        view(3)
        xlabel('tSNE1');
        ylabel('tSNE2');
        zlabel('tSNE3');
        title(genes_to_show{i});
        Ax = [Ax, ax];

        % tSNE colored by RGB
        ax = subplot(223); hold on;
        for j = 1:size(pos,1)
            plot3(tsne_3d(j,1), tsne_3d(j,2), tsne_3d(j,3), '.',...
                'color', tsne_rgb(j,:));
        end
        view(3)
        xlabel('tSNE1');
        ylabel('tSNE2');
        zlabel('tSNE3');
        title('tSNE in RGB');
        Ax = [Ax, ax];

        % spatial map colored by expression
        subplot(222); hold on;
        for j = 1:size(pos,1)
            [vy, vx] = dot2poly(pos(j,2), pos(j,1), hexbin_size, 6);
            patch(vx, vy, col(counts(j,idx)+1,:));
        end
        axis image
        axis off
        title(genes_to_show{i});
        ch = colorbar(gca, 'Colormap', col);
        ch.Ticks = [0 1];
        ch.TickLabels = [0 max(counts(:,idx))];
        ch.Label.String = 'count in a bin';

        % spatial map colored by tSNE RGB
        subplot(224); hold on;
        for j = 1:size(pos,1)
            [vy, vx] = dot2poly(pos(j,2), pos(j,1), hexbin_size, 6);
            patch(vx, vy, tsne_rgb(j,:));
        end
        axis image
        axis off
        title('RGB in spatial');

        % synchronize two tSNE plots
        linkprop(Ax, {'view', 'xlim', 'ylim', 'zlim'});

        % wait for keyboard input to continue
        pause()
    end
end

