function plot_all_onblack(infile)
% quickly plot all reads on a black background, with axis
% Xiaoyan, 2017


% input and remove NNNN
[name, pos] = getinsitudata(infile);
[name, pos] = removereads(name, 'NNNN', pos);

plotall(name, pos);

axis on
set(gca, 'color', 'k');

end
