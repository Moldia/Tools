[names,pos] = getinsitudata_f('E:\Whole_organoid_pseudoAnchor_v5\Decoding\QT_0.7_0_details.csv');
[name_uni,~,idx_name] = unique(names);
idx_NNNN = find(strcmp('NNNN',name_uni));

names = names(idx_name~=idx_NNNN);
pos = pos(idx_name~=idx_NNNN,:);
[name_uni,~,idx_name] = unique(names);
name_count = hist(idx_name,1:length(name_uni));

NN = rangesearch(pos,pos,5);
NN_max = max(cellfun(@(v) length(v),NN));
NN = cellfun(@(v) [v, nan(1,NN_max-length(v))],NN,'uni',0);
NN = cell2mat(NN);
NN = NN(:,2:end);
NNmat = zeros(length(name_uni),length(name_uni),'uint8');
for i = 1:length(name_uni)
    temp = idx_name == i;
    temp = NN(temp,:);
    temp = temp(:);
    temp(isnan(temp)) = [];
    temp = idx_name(temp);
    NNmat(i,:) = hist(temp,1:length(name_uni));
end


%% simulations
NNmat_rand = zeros(length(name_uni),length(name_uni),5,'uint8');

for r = 1:5
    r
    randidx = randperm(length(names));
    randnames = names(randidx);
    
    [name_uni,~,idx_name] = unique(randnames);
    name_count = hist(idx_name,1:length(name_uni));
    
    NN = rangesearch(pos,pos,5);
    NN_max = max(cellfun(@(v) length(v),NN));
    NN = cellfun(@(v) [v, nan(1,NN_max-length(v))],NN,'uni',0);
    NN = cell2mat(NN);
    NN = NN(:,2:end);
    for i = 1:length(name_uni)
        temp = idx_name == i;
        temp = NN(temp,:);
        temp = temp(:);
        temp(isnan(temp)) = [];
        temp = idx_name(temp);
        NNmat_rand(i,:,r) = hist(temp,1:length(name_uni));
    end
end

NNmat_rand = mean(NNmat_rand,3);
%%
figure,
Ax = [];

ax = subplot(1,2,1);
bh = bar3(NNmat);
for i = 1:length(bh)
    bh(i).CData = bh(i).ZData;
end
view(2)
set(gca,'xtick',1:length(name_uni),'xticklabel',name_uni,...
    'ytick',1:length(name_uni),'yticklabel',name_uni,...
    'xticklabelrotation',90)
colorbar
daspect([1 1 1])
Ax = [Ax, ax];

ax = subplot(1,2,2);
bh = bar3(NNmat_rand);
for i = 1:length(bh)
    bh(i).CData = bh(i).ZData;
end
view(2)
set(gca,'xtick',1:length(name_uni),'xticklabel',name_uni,...
    'ytick',1:length(name_uni),'yticklabel',name_uni,...
    'xticklabelrotation',90)
colorbar
daspect([1 1 1])
Ax = [Ax, ax];

% ax = subplot(2,2,2);
% bh = bar3(double(NNmat)./repmat(name_count',1,length(name_uni)));
% for i = 1:length(bh)
%     bh(i).CData = bh(i).ZData;
% end
% view(2)
% set(gca,'xtick',1:length(name_uni),'xticklabel',name_uni,...
%     'ytick',1:length(name_uni),'yticklabel',name_uni,...
%     'xticklabelrotation',90)
% colorbar
% daspect([1 1 1])
% Ax = [Ax, ax];
% 
% ax = subplot(2,2,3);
% bh = bar3(sqrt(double(NNmat)./repmat(name_count',1,length(name_uni))./repmat(name_count,length(name_uni),1)));
% for i = 1:length(bh)
%     bh(i).CData = bh(i).ZData;
% end
% view(2)
% set(gca,'xtick',1:length(name_uni),'xticklabel',name_uni,...
%     'ytick',1:length(name_uni),'yticklabel',name_uni,...
%     'xticklabelrotation',90)
% colorbar
% daspect([1 1 1])
% Ax = [Ax, ax];
% 
% ax = subplot(2,2,4);
% bh = bar3(double(NNmat)./(repmat(name_count',1,length(name_uni))+repmat(name_count,length(name_uni),1)));
% for i = 1:length(bh)
%     bh(i).CData = bh(i).ZData;
% end
% view(2)
% set(gca,'xtick',1:length(name_uni),'xticklabel',name_uni,...
%     'ytick',1:length(name_uni),'yticklabel',name_uni,...
%     'xticklabelrotation',90)
% colorbar
% daspect([1 1 1])
% Ax = [Ax, ax];


linkaxes(Ax,'xy');
