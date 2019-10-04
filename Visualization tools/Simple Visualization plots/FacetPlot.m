% plot each transcript separately
% Xiaoyan, 2017

close all; drawnow;

%% modify here
decoded_file = 'E:\PROOOJECTS\test_dataset\QT_0.35_0.004_details.csv';
output_folder = 'E:\PROOOJECTS\test_dataset\FacePlot'; 

%% do not modify

% load and remove NNNN
[name, pos] = getinsitudata(decoded_file);
[name, pos] = removereads(name, 'NNNN', pos);

% unique transcripts
[uNames, ~, iName] = unique(name);
cName = hist(iName, 1:length(uNames));

% plot
f = 0;
for i = 1:numel(uNames)
    if mod(i-1, 20)==0
        figure;
        f = f+1;
        set(gcf, 'units', 'normalized',...
            'position', [.05+.02*f .1-.02*f .8 .8],...
            'visible','off');
        s = 0;
    end
    
    s = s+1;
    subplot(4,5,s);
    plotonblank;
    plot(pos(iName==i,1), pos(iName==i,2), '.');
    title([uNames{i} ' (' num2str(cName(i)) ')']);
end

% save output
try
    mkdir(output_folder);
end

disp('saving images..');
while f>=1
    figure(f);
    set(gcf, 'visible', 'on');
    saveas(gcf, fullfile(output_folder, ['FacetPlot_' num2str(f) '.png']), 'png');
    f = f-1;
end
