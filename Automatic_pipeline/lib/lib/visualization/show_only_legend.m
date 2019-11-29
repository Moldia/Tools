function show_only_legend(ah, legendnames, symbols, fontsize)
% show_only_legend(ah, legendnames, symbols)
% create an invisible dummy figure in subplot to show only legend
% Xiaoyan, 2019


% randomPoints = rand(length(legendnames),2);

axes(ah);
hold on;
for i = 1:numel(legendnames)
    ph = plot(nan, nan, symbols{i});
%     set(ph, 'visible', 'off');
end
if nargin > 3
    h = legend(legendnames, 'color', [.6 .6 .6], 'fontsize', fontsize);
else
    h = legend(legendnames, 'color', [.6 .6 .6], 'fontsize', 5);    
end
set(h, 'location', 'west');
set(ah, 'visible', 'off');

end
