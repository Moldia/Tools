%% compare relative closeness of A and B to any other transcript
%  Xiaoyan 2014-11-17

%% import data
data = importdata('E:\Whole_organoid_pseudoAnchor_v5\Decoding\QT_0.7_0_details.csv',',',1);
pos = data.data(:,1:2);
name = data.textdata(2:end,2);
clear data

%% unique transcripts
[name_uni,~,re_idx] = unique(name);
idx_A = find(strcmp('Homo sapiens microtubule associated protein 2 (MAP2)',name_uni));
idx_B = find(strcmp('Homo sapiens tubulin beta 3 class III (TUBB3)',name_uni));
idx_NNNN = find(strcmp('NNNN',name_uni));
idx_rest = 1:length(name_uni);
idx_rest([idx_A,idx_B,idx_NNNN]) = [];

[a,b]=hist(re_idx,unique(re_idx));
figure;
bar(b,a);
set(gca,'XTick',1:length(name_uni),...
    'XTickLabel',name_uni,'XTickLabelRotation',90,...
    'XLim',[0 length(name_uni)+1]);

%% plot receptor A and B transcripts
radius = 150;
alph = linspace(0,2*pi,21);
circ_x = radius*cos(alph);
circ_y = radius*sin(alph);

pos_A = pos(re_idx==idx_A,:);
pos_B = pos(re_idx==idx_B,:);

% poly_Ax = repmat(pos_A(:,1),1,21) + repmat(circ_x,size(pos_A,1),1);
% poly_Ay = repmat(pos_A(:,2),1,21) + repmat(circ_y,size(pos_A,1),1);
% poly_Bx = repmat(pos_B(:,1),1,21) + repmat(circ_x,size(pos_B,1),1);
% poly_By = repmat(pos_B(:,2),1,21) + repmat(circ_y,size(pos_B,1),1);
% figure; hold on;
% axis image
% set(gca,'YDir','reverse');
% patch(poly_Ax',poly_Ay','r','facealpha',.3);
% patch(poly_Bx',poly_By','b','facealpha',.3);

%% find neighbors of transcripts except for receptor A and B
pos_pool = [pos_A,ones(size(pos_A,1),1);pos_B,repmat(2,size(pos_B,1),1)];
% poly_poolx = repmat(pos_pool(:,1),1,21) + repmat(circ_x,size(pos_pool,1),1);
% poly_pooly = repmat(pos_pool(:,2),1,21) + repmat(circ_y,size(pos_pool,1),1);

NN_hist = [];
NN_idx = [];
D = [];
for i = 1:length(name_uni)
    i
    if i == idx_A || i == idx_B || i == idx_NNNN
    else
        pos_q = pos(re_idx==i,:);
%         poly_qx = repmat(pos_q(:,1),1,21) + repmat(circ_x,size(pos_q,1),1);
%         poly_qy = repmat(pos_q(:,2),1,21) + repmat(circ_y,size(pos_q,1),1);
%         patch(poly_qx',poly_qy','c','facealpha',.3);
        
        [IDX,d] = knnsearch(pos_pool(:,1:2),pos_q,...
            'NSmethod','kdtree','Distance','euclidean');
%         patch(poly_poolx(I,:)',poly_pooly(I,:)','y','facealpha',.3);
%         hold on;
%         patch(poly_qx',poly_qy','r','facealpha',.3);
        
        NN_hist = [NN_hist;hist(pos_pool(IDX,3),1:2)];   
        NN_idx = [NN_idx;{IDX}];
        D = [D;mean(d(pos_pool(IDX,3)==1)),mean(d(pos_pool(IDX,3)==2))];
    end
end

NN_hist = NN_hist./repmat(sum(NN_hist,2),1,2);
random_AB = size(pos_A,1)/size(pos_pool,1);

%% bar plot
% figure;bar(NN_hist,'stacked');
% set(gca,'XTick',1:length(name_rest),...
%     'XTickLabel',name_uni(name_rest),'XTickLabelRotation',90,...
%     'XLim',[0 length(name_rest)+1]);
% 
% hold on;
% plot([0 length(name_rest)+1],...
%     [random_AB,random_AB],...
%     'r','linewidth',2);
% plot([0 length(name_rest)+1],...
%     [random_AB,random_AB]-.05,...
%     'r','linewidth',1,'linestyle','--');
% plot([0 length(name_rest)+1],...
%     [random_AB,random_AB]+.05,...
%     'r','linewidth',1,'linestyle','--');
% plot([0 length(name_rest)+1],...
%     [random_AB,random_AB]-.1,...
%     'r','linewidth',.5,'linestyle','-.');
% plot([0 length(name_uni)-1],...
%     [random_AB,random_AB]+.1,...
%     'r','linewidth',.5,'linestyle','-.');


%% plot specific transcript - PDGFR
figure; hold on;
axis image
set(gca,'YDir','reverse');

name_test = find(strcmp('Homo sapiens doublecortin (DCX)',name_uni));
pos_test = pos(re_idx==name_test,:);
poly_testx = repmat(pos_test(:,1),1,21) + repmat(circ_x,size(pos_test,1),1);
poly_testy = repmat(pos_test(:,2),1,21) + repmat(circ_y,size(pos_test,1),1);
patch(poly_testx',poly_testy','g','facealpha',.9,'linestyle','none');

list = 1:length(name_uni);
list([idx_A,idx_B,idx_NNNN]) = [];
idx_test = find(list==name_test);
pos_NN = NN_idx{idx_test};
pos_NN = pos_pool(pos_NN,:);

pos_NN_A = pos_NN(pos_NN(:,3)==1,1:2);
pos_NN_B = pos_NN(pos_NN(:,3)==2,1:2);
poly_NN_Ax = repmat(pos_NN_A(:,1),1,21) + repmat(circ_x,size(pos_NN_A,1),1);
poly_NN_Ay = repmat(pos_NN_A(:,2),1,21) + repmat(circ_y,size(pos_NN_A,1),1);
poly_NN_Bx = repmat(pos_NN_B(:,1),1,21) + repmat(circ_x,size(pos_NN_B,1),1);
poly_NN_By = repmat(pos_NN_B(:,2),1,21) + repmat(circ_y,size(pos_NN_B,1),1);
patch(poly_NN_Ax',poly_NN_Ay','r','facealpha',.3,'linestyle','none');
patch(poly_NN_Bx',poly_NN_By','b','facealpha',.3,'linestyle','none');

%% neighbors of PDGFRA and PDGFRB separately
% NN_idx_A = [];
D_A = [];
for i = 1:length(name_uni)
    if i == idx_A || i == idx_B || i == idx_NNNN
    else
        pos_q = pos(re_idx==i,:);
        
        [IDX,d] = knnsearch(pos_A(:,1:2),pos_q,...
            'NSmethod','kdtree','Distance','euclidean');
        
%         NN_idx_A = [NN_idx_A;length(I)];   
        D_A = [D_A;name_uni(i),mean(d),std(d)];
    end
end

% NN_idx_B = [];
D_B = [];
for i = 1:length(name_uni)
    if i == idx_A || i == idx_B || i == idx_NNNN
    else
        pos_q = pos(re_idx==i,:);
        
        [IDX,d] = knnsearch(pos_B(:,1:2),pos_q,...
            'NSmethod','kdtree','Distance','euclidean');
        
%         NN_idx_B = [NN_idx_B;length(I)];   
        D_B = [D_B;name_uni(i),mean(d),std(d)];
    end
end

% figure; 
% subplot(3,1,[1,2]);hold on;
% for i = 1:length(name_rest)
%     plot([i,i],[D_A{i,2},D_B{i,2}],'c');
% end
% plot(1:length(name_rest),cell2mat(D_A(:,2)),'.');
% plot(1:length(name_rest),cell2mat(D_B(:,2)),'.');
% 
% subplot(3,1,3);
% bar(cell2mat(D_B(:,2))-cell2mat(D_A(:,2)))
% set(gca,'XTick',1:length(name_rest),...
%     'XTickLabel',name_uni(name_rest),'XTickLabelRotation',90,...
%     'XLim',[0 length(name_rest)+1]);

%% A/B ratio of neighboring transcripts
NN_hist_rate = [];
tic
for r = 100:100:300
    hist_rate_r = [];
    for i = 1:length(name_uni)
        i
        if i == idx_A || i == idx_B || i == idx_NNNN
        else
            pos_q = pos(re_idx==i,:);
            
            IDX = rangesearch(pos_pool(:,1:2),pos_q,r,...
                'NSmethod','kdtree','Distance','euclidean');
            
            num = cellfun(@(v) size(v),IDX,'uniformoutput',false);
            num = cell2mat(num);
            num = num(:,2);
            [num_uni,~,num_re] = unique(num);
            
            I = [];
            for j = 1:length(num_uni)
                temp = cell2mat(IDX(num_re==j));
                I = [I;temp(:)];
            end
            
            hist_rate = hist(pos_pool(I,3),1:2);
            hist_rate_r = [hist_rate_r;{hist_rate}];
        end
    end
    NN_hist_rate = [NN_hist_rate,[{num2str(r)};hist_rate_r]];
end
toc

figure;
set(gca,'XTick',1:length(idx_rest),...
    'XTickLabel',name_uni(idx_rest),'XTickLabelRotation',90,...
    'XLim',[0 length(idx_rest)+1],'YLim',[0 1]);
h2 = axes('Position',get(gca,'Position'));
set(h2,'XAxisLocation','top','YTick',[],'YLim',[0 1],...
    'XTick',1:length(idx_rest),...
    'XTickLabel',a(idx_rest),'XTickLabelRotation',90)
hold on;
for r = 1:size(NN_hist_rate,2)
    temp = cell2mat(NN_hist_rate(2:end,r));
    temp = temp./repmat(sum(temp,2),1,2);
    plot(1:length(idx_rest),temp(:,1),'.','markersize',15)
end
plot([0 length(idx_rest)+1],...
    [random_AB,random_AB],...
    'r','linestyle','--','linewidth',.5);
legend(NN_hist_rate(1,:));
title('K\_seq');

save('Neighbor.mat','a','random_AB','NN_idx','name_uni','NN_hist_rate');

