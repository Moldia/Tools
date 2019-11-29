function [uniName, probenum] = find_oligo_genes(fileUsedProbes, allProbes)
% [uniName, probenum] = find_oligo_genes(fileUsedProbes, allProbes)
% input: file containg list of oligo ID, corresponding search library
% output: targets and number of probes
% Xiaoyan, 2017


queryProbes = importdata(fileUsedProbes);
hits = cellfun(@(v) strcmp(v, allProbes(:,1)), queryProbes, 'uni', 0);
hits = cellfun(@find, hits);    

genes = allProbes(hits,2);
[uniName, ~, idxName] = unique(genes);

probenum = hist(idxName, 1:length(uniName));

end