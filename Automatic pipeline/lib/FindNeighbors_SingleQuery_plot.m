%% written fro Oncotrak project
%  find stroma/ tumor neighbors of Igf2


%% list of files to read
files = ls;
files = files(find(strcmp('QT',cellstr(files(3:end,1:2))))+2,:);
files = cellstr(files);

%% find the number of Igf2 found close to itself, to stroma, to tumor within a certain distance
distance = 15;  % unit: pixel

NNmat = zeros(length(files),3);
NNcoord = cell(length(files),3);
Igf2_TotalCount = zeros(length(files),1);

for i = 1:length(files)
    [names,pos] = getinsitudata_f(files{i});
    [name_uni,~,idx_name] = unique(names);
    
    % remove NNNN
    idx_NNNN = find(strcmp('NNNN',name_uni));
    names = names(idx_name~=idx_NNNN);
    pos = pos(idx_name~=idx_NNNN,:);
    [name_uni,~,idx_name] = unique(names);

    % count Igf2 transcripts
    name_count = hist(idx_name,1:length(name_uni));
    idx_Igf2 = find(strcmp('Igf2',name_uni));
    Igf2_TotalCount(i) = name_count(idx_Igf2);

    % range search
    NN = rangesearch(pos,pos,distance);
    NN_max = max(cellfun(@(v) length(v),NN));
    NN = cellfun(@(v) [v, nan(1,NN_max-length(v))],NN,'uni',0);
    NN = cell2mat(NN);
    NN = NN(:,2:end);
 
    % Igf2 as query
    temp = idx_name == idx_Igf2;
    coord_q = pos(temp,:);

    % connect query and hits, find query indeces and coordinates
    temp = NN(temp,:);
    temp_q = double(temp>0).*repmat((1:length(temp))',1,NN_max-1);
    temp_q = temp_q(:);
    temp_q(temp_q==0) = [];
    coord_q = coord_q(temp_q,:);
 
    % remove Igf2 singlets, find hit indeces and coordinates
    temp = temp(:);
    temp(isnan(temp)) = [];
    coord_h = pos(temp,:);
    
    % map transcripts, count number of hits
    temp = idx_name(temp);
    NNmat(i,:) = hist(temp,1:length(name_uni))';
    
    % store query-hit coordinates pairs for each hit transcript species
    for j = 1:length(name_uni)
        coord_pair = [coord_q(temp==j,:),coord_h(temp==j,:)];
        NNcoord{i,j} = coord_pair;
    end
end

%% plot
mkdir('plot');
sym = {'ro','go','bo'};
images = {'p150_HE_s2_sc100.png','p150_HE_s4m1_sc100.png','p150_HE_s4m2_sc100.png','p150_HE_s8_sc90.png','p150_HE_s6_sc100.png'};
scales = [1 1 1 0.9 1];
for i = 1`:length(files)
    clf;
    imshow(images{i});
    hold on;
    plot(NNcoord{i,2}(:,1)*scales(i),NNcoord{i,2}(:,2)*scales(i),'k+');
    plot(NNcoord{i,2}(:,3)*scales(i),NNcoord{i,2}(:,4)*scales(i),'bo');
    plot(NNcoord{i,3}(:,3)*scales(i),NNcoord{i,3}(:,4)*scales(i),'ro');
    plot(NNcoord{i,3}(:,1)*scales(i),NNcoord{i,3}(:,2)*scales(i),'k+');
    for j = 1:length(NNcoord{i,2})
        plot([NNcoord{i,2}(j,1),NNcoord{i,2}(j,3)]*scales(i),[NNcoord{i,2}(j,2),NNcoord{i,2}(j,4)]*scales(i),'k-');
    end
    for j = 1:length(NNcoord{i,3})
        plot([NNcoord{i,3}(j,1),NNcoord{i,3}(j,3)]*scales(i),[NNcoord{i,3}(j,2),NNcoord{i,3}(j,4)]*scales(i),'k-');
    end
%     axis image;
%     set(gca,'xlim',[0 max(pos(:,1))],'ylim',[0 max(pos(:,2))],'YDir','reverse')
%     axis off
    legend(name_uni,'location','northeastoutside');
    drawnow;
%     print(gcf,['plot\' files{i},'.png'],'-dpng');
    saveas(gcf,['plot\' files{i} '.fig']);
end

%% format into table, and write
table_NN = array2table(NNmat,'VariableNames',name_uni');
table_NN = [table_NN,table(Igf2_TotalCount)];
table_NN.Properties.RowNames = files;
writetable(table_NN,'number_of_Igf2_found_close_to.csv','WriteRowNames',True);

figure,bar(bsxfun(@rdivide,NNmat,Igf2_TotalCount)');
set(gca,'xticklabel',name_uni);
legend(files);