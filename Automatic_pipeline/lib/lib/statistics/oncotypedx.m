function [RiskGroup, Counts] = oncotypedx(name)
%% calculate OncotypeDX score
%  Xiaoyan, 2016-9-15

%% import
[name_uni, ~, idx_name] = unique(name);

%% gene panel for OncotypeDX
Genes.Ref = {'ACTB','GAPDH','RPLP0','GUS1','TFRC'};
Genes.Prolif = {'Ki-67','STK15','BIRC5','CCNB1','MYBL2'};   % survivin aka BIRC5
Genes.Invasion = {'MMP11','CTSL2'};
Genes.GRB7 = {'GRB7','HER2'};
Genes.ER = {'ER','PR','BCL2','SCUBE2'};
Genes.Other = {'GSTM1','CD68','BAG1'};

%% factor for scoring
Factor.GRB7 = [.9, .1];
Factor.ER = [.8, 1.2, 1, 1]/4;
Factor.Prolif = [1, 1, 1, 1, 1]/5;
Factor.Invasion = [1, 1]/2;
Factor.Other = [-.08, .05, -.07];
RSUCoeff = [.47, -.34, 1.04, .1, 1];

%% find genes
Groups = fieldnames(Genes);
IdxGenes = struct();

for i = 1:length(Groups)
    tempname = Genes.(Groups{i});
    tempidx = [];
    for j = 1:length(tempname)
        try
            tempidx(j) = find(strcmp(tempname{j}, name_uni));
        catch ME
            if strcmp(ME.identifier, 'MATLAB:badRectangle')		% catch error: dimension mismatch
                tempidx(j) = 0;		% no read found for this gene
            end
        end
    end
    IdxGenes.(Groups{i}) = tempidx;
end

%% normalization
Counts = struct();
CountsAll = 0;
for i = 1:length(Groups)
    tempidx = IdxGenes.(Groups{i});
    tempcount = [];
    for j = 1:length(tempidx)
        tempcount = [tempcount, nnz(idx_name==tempidx(j))];
        CountsAll = CountsAll + nnz(idx_name==tempidx(j));
    end
    Counts.(Groups{i}) = tempcount;
end

AverageRef = sum(Counts.Ref)/5;

%% classification
if AverageRef>0 && CountsAll>20

    CountsNorm = structfun(@(v) abs(log2((v+1)/AverageRef)), Counts, 'uni', 0);

    %% normalized values range from 0 to 15
    CountsNorm = structfun(@(v) max(0,v), CountsNorm, 'uni', 0);
    CountsNorm = structfun(@(v) min(15,v), CountsNorm, 'uni', 0);

    %% group scores
    GroupScore = zeros(6,1);
    for i = 2:length(Groups)
        GroupScore(i) = sum(Factor.(Groups{i}).*CountsNorm.(Groups{i}));
    end

%     if GroupScore(strcmp(Groups,{'GRB7'})) < 8
%         GroupScore(strcmp(Groups,{'GRB7'})) = 8;
%     end
%     
%     if GroupScore(strcmp(Groups,{'Prolif'})) < 6.5
%         GroupScore(strcmp(Groups,{'Prolif'})) = 6.5;
%     end
    
    %% unscaled recurrence score
    RSU = sum(RSUCoeff'.*GroupScore(2:end));

    %% scaling
    if RSU < 0
        RS = 0;
    elseif RSU <= 100
        RS = 20*(RSU-6.7);
    else
        RS = 100;
    end

    %% classification
    if RS < 18
        RiskGroup = 'low';
    elseif RS < 31
        RiskGroup = 'intermediate';
    else
        RiskGroup = 'high';
    end
else
    RiskGroup = 'NA';
end

end