% plot selected subgroup of transcritps
% Xiaoyan, 2018


%% modify here
decoded_file = 'E:\PROOOJECTS\test_dataset\QT_0.35_0.004_details.csv';
image = 'E:\PROOOJECTS\test_dataset\860502_1_align.png';
scale = .2; % image scale
genes_to_show = {'VIM', 'ACTB'};
figure_title = 'whatever markers';

%% do not modify

% load
[name, pos] = getinsitudata(decoded_file);

% plot
figure;
plotall(name, pos, image, scale);
update_legend(gca, genes_to_show);

