% test the null hypothesis if two groups of data points are drawn from the same continuous distribution 
% two-sample Kolmogorov-Smirnov test
% Xiaoyan, 2017


input_file = 'E:\Organoid 1\Organoid 1S\QT_0.45_1e-07_details.csv';
% genes with read lower than this threshold will not be used in analysis
count_threshold = 30; 
% size of the bin
bin_size = 200;


%%
[name, pos] = getinsitudata(input_file);
[uNames, ~, idxName] = unique(name);
cNames = hist(idxName, 1:max(idxName));

% remove NNNN and any reads<count_threshold
[name, pos] = removereads(name, [uNames(cNames<count_threshold); {'NNNN'}], pos);
[uNames, ~, idxName] = unique(name);

% bin
plot_all_onblack(input_file);
show_square_lines(gcf, bin_size, bin_size, 'w');
pos_bin = ceil(pos/bin_size);
pos_bin = max(pos_bin(:,1))*(pos_bin(:,1)-1) + pos_bin(:,2);


ksstat = zeros(length(uNames), 2);
for i = 1:length(uNames)
    for j = i:length(uNames)
        [~, p, ks2stat] = kstest2(pos_bin(idxName==i), pos_bin(idxName==j));
        ksstat(i,j,1) = ks2stat;
        ksstat(i,j,2) = p;
        ksstat(j,i,1) = ks2stat;
        ksstat(j,i,2) = p;
    end
end
% figure; imagesc(ksstat(:,:,2));
% set(gca, 'xtick', 1:length(uniName), 'xticklabel', uniName,...
%     'ytick', 1:length(uniName), 'yticklabel', uniName,...
%     'xticklabelrotation', 90);


% sorting
L = linkage(ksstat(:,:,2));
order = dendroperm(L);

figure; imagesc(ksstat(order,order,2)); grid on;
set(gca, 'xtick', 1:length(uNames), 'xticklabel', uNames(order),...
    'ytick', 1:length(uNames), 'yticklabel', uNames(order),...
    'xticklabelrotation', 90);

c = colorbar;
c.Label.String = 'Kolmogorov-Smirnov test p-value';
c.Label.FontSize = 12;
