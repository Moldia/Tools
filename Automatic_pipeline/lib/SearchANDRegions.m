% search regions where all specified genes co-occur
% Xiaoyan, 2017

decoded_file = 'E:\Whole_organoid_pseudoAnchor_v5\Decoding\QT_0.7_0_details.csv';
image = 'K:\Organoid_method_development\Organoids_erik\Organoid3_Output\3\Sample 3_Base 3_resize20_c1.jpg';
scale = 0.2;

%%
[name, pos] = getinsitudata(decoded_file);
[name, pos] = removereads(name, 'NNNN', pos);
figure;
plotall(name, pos, image, scale)


 %cooccur = search_reads_cooccur(name, pos, 100);
% cooccur = search_reads_cooccur(name, pos, 100, {'Drd3', 'Prrx1', 'Ctdnep1'});
cooccur = search_reads_combs(name, pos, 100,2,{'Homo sapiens tubulin beta 3 class III (TUBB3)', 'Homo sapiens doublecortin (DCX)'});

