% search regions where all specified genes co-occur
% Xiaoyan, 2017

decoded_file = 'Decoding\QT_0.45_0_details.csv';
image = 'slideA_alignedHE.jpg';
scale = 1;

%%
[name, pos] = getinsitudata(decoded_file);
[name, pos] = removereads(name, 'NNNN', pos);
figure;
plotall(name, pos, image, scale)


% cooccur = search_reads_cooccur(name, pos, 100);
% cooccur = search_reads_cooccur(name, pos, 100, {'Drd3', 'Prrx1', 'Ctdnep1'});
cooccur = search_reads_combs(name, pos, 100, 3, {'Drd3', 'Prrx1'});

