%% written by Christoffer M. Langseth and Markus M. Hilscher, 2019
% collection of scripts for compartive analysis of two image analysis
% pipelines
%% read files
nilsson = readtable('/Volumes/Christoffer/SCRIPTSnilssonStarfishComp/nilsson.csv');
%nilsson_sorted = sortrows(nilsson,3,'ascend'); % ordering the coordinates 
starfish = readtable('/Volumes/Christoffer/SCRIPTSnilssonStarfishComp/starfish.csv'); 
%starfish_sorted = sortrows(starfish,15,'ascend'); % ordering the coordinates
%% defining x, y and the grouping

% before normalization and conversion
nilsson_x_sorted = nilsson_sorted.PosX;
nilsson_y_sorted = nilsson_sorted.PosY;
nilsson_g_sorted = nilsson_sorted.Gene;

starfish_x_sorted = starfish_sorted.xc;
starfish_y_sorted = starfish_sorted.yc;
starfish_g_sorted = starfish_sorted.Gene;

% global coordinates (um) to coordinates in pixel space
starfish_x_sorted_normalized = (starfish_sorted.xc/0.1625);
starfish_y_sorted_normalized = (starfish_sorted.yc/0.1625);

% correcting for the shift from the automated alignment
starfish_x_sorted_normalized_shift = (starfish_sorted.xc/0.1625)+18.7911;
starfish_y_sorted_normalized_shift = (starfish_sorted.yc/0.1625)-3.2584;

csvwrite('starfish_sorted.csv', starfish_x_sorted_normalized_shift);
csvwrite('starfish_sorted.csv', starfish_y_sorted_normalized_shift)


%% side by side compairson: plotting the coordinates 

% before normalization and conversion;
figure();
title("Gene maps created by the Nilsson and the Starfish pipeline, before normalization and conversion")
subplot(1,2,1)
gscatter(nilsson_x_sorted, nilsson_y_sorted, nilsson_g_sorted, '', '', 2, 'on','x', 'y');

subplot(1,2,2)
gscatter(starfish_x_sorted, starfish_y_sorted, starfish_g_sorted, '', '', 2,'on', 'x', 'y');

% global coordinates (um) to coordinates in pixel space
figure();
title("Gene maps created by the Nilsson and the Satrfish pipeline, pixel")
subplot(1,2,1)
gscatter(nilsson_x_sorted, nilsson_y_sorted, nilsson_g_sorted, '', '', 2, 'on','x', 'y');

subplot(1,2,2)
gscatter(starfish_x_sorted_normalized, starfish_y_sorted_normalized, starfish_g_sorted, '', '', 2,'on', 'x', 'y');

% correcting for the shift from the automated alignment
figure();
title("Gene maps created by the Nilsson and the Satrfish pipeline, ")
subplot(1,2,1)
gscatter(nilsson_x_sorted, nilsson_y_sorted, nilsson_g_sorted, '', '', 2, 'on','x', 'y');

subplot(1,2,2)
gscatter(starfish_x_sorted_normalized_shift, starfish_y_sorted_normalized_shift, starfish_g_sorted, '', '', 2,'on', 'x', 'y');


%% overlaying 

% before conversion
figure() 
gscatter(nilsson_x_sorted, nilsson_y_sorted, nilsson_g_sorted, 'g', '*', 3);
hold on
gscatter(starfish_x_sorted, starfish_y_sorted, starfish_g_sorted, 'r', 'o', 3);
hold on
title('Overlayed non-normalized values')
xlabel('x')
ylabel('y')
ax = gca;
ax.FontSize = 16;

% global coordinates to pixel space
figure() 
gscatter(nilsson_x_sorted, nilsson_y_sorted, nilsson_g_sorted, 'g', '*', 3);
hold on
gscatter(starfish_x_sorted_normalized, starfish_y_sorted_normalized, starfish_g_sorted, 'r', 'o', 3);
hold on
title('Overlayed global coordinates to pixel space-converted values')
xlabel('x')
ylabel('y')
ax = gca;
ax.FontSize = 16;

% correcting for shift
figure() 
gscatter(nilsson_x_sorted, nilsson_y_sorted, nilsson_g_sorted, 'g', '*', 3);
hold on
gscatter(starfish_x_sorted_normalized_shift, starfish_y_sorted_normalized_shift, starfish_g_sorted, 'r', 'o', 3);
hold on
title('Overlayed shifted values')
xlabel('x')
ylabel('y')
ax = gca;
ax.FontSize = 16;

%% rounding
% global coordinates (um) to pixel space-converted  
nilsson_x_sorted_rounded = round(nilsson_x_sorted, -2);
nilsson_y_sorted_rounded = round(nilsson_y_sorted, -2);

starfish_x_sorted_normalized_rounded = round(starfish_x_sorted_normalized, -2);
starfish_y_sorted_normalized_rounded = round(starfish_y_sorted_normalized, -2);

%global coordinates (um) to pixel space-converted and shifted
nilsson_x_sorted_rounded = round(nilsson_x_sorted, -2);
nilsson_y_sorted_rounded = round(nilsson_y_sorted, -2);

starfish_x_sorted_normalized_shift_rounded = round(starfish_x_sorted_normalized_shift, -2);
starfish_y_sorted_normalized_shift_rounded = round(starfish_y_sorted_normalized_shift, -2);

% % testing out other degrees of rounding 
% nilsson_x_sorted_rounded = round(nilsson_x_sorted, -1);
% nilsson_y_sorted_rounded = round(nilsson_y_sorted, -1);
% 
% starfish_x_sorted_rounded = round(starfish_x_sorted, -1);
% starfish_y_sorted_rounded = round(starfish_y_sorted, -1);
% 
% nilsson_x_sorted_rounded = round(nilsson_x_sorted, -3);
% nilsson_y_sorted_rounded = round(nilsson_y_sorted, -3);
% 
% starfish_x_sorted_rounded = round(starfish_x_sorted, -3);
% starfish_y_sorted_rounded = round(starfish_y_sorted, -3);

% nilsson_x_sorted_rounded = round(nilsson_x_sorted, -4);
% nilsson_y_sorted_rounded = round(nilsson_y_sorted, -4);
% 
% starfish_x_sorted_rounded = round(starfish_x_sorted, -4);
% starfish_y_sorted_rounded = round(starfish_y_sorted, -4);

%% 
% % rounded values
% N_2=[nilsson_x_sorted_normalized_rounded, nilsson_y_sorted_normalized_rounded];
% S_2=[starfish_x_sorted_normalized_rounded, starfish_y_sorted_normalized_rounded];
% 
% [C_2,iN_2,iS_2] = intersect(N_2,S_2,'rows');
% 
% % pixel normalized and automated image registration normalized
% N_2=[nilsson_x_sorted, nilsson_y_sorted];
% S_2=[starfish_x_sorted, starfish_y_sorted];
% 
% [C_2,iN_2,iS_2] = intersect(N_2,S_2,'rows');

% % only global coordinates (um) to pixel space-converted
% N_2=[nilsson_x_sorted, nilsson_y_sorted];
% S_2=[starfish_x_sorted_normalized, starfish_x_sorted_normalized];
% [C_2,iN_2,iS_2] = intersect(N_2,S_2,'rows');
% 
% % only global coordinates (um) to pixel space-converted and rounded
% N_2=[nilsson_x_sorted_rounded, nilsson_y_sorted_rounded];
% S_2=[starfish_x_sorted_normalized_rounded, starfish_y_sorted_normalized_rounded];
% [C_2,iN_2,iS_2] = intersect(N_2,S_2,'rows');
% 
% % only global coordinates (um) to pixel space-converted and shifted
% N_2=[nilsson_x_sorted_rounded, nilsson_y_sorted_rounded];
% S_2=[starfish_x_sorted_normalized_shift, starfish_y_sorted_normalized_shift];
% [C_2,iN_2,iS_2] = intersect(N_2,S_2,'rows');

% global coordinates (um) to pixel space-converted, shifted rounded
N_2=[nilsson_x_sorted_rounded, nilsson_y_sorted_rounded];
S_2=[starfish_x_sorted_normalized_shift_rounded, starfish_y_sorted_normalized_shift_rounded];
[C_2,iN_2,iS_2] = intersect(N_2,S_2,'rows');
%%
nilsson_genes = nilsson_sorted(iN_2,3);
starfish_genes = starfish_sorted(iS_2,15);

%nilsson_genes = starfish(iN_2,15);
%starfish_genes = nilsson(iN_2,3);

nilsson_genes = table2array(nilsson_genes);
starfish_genes = table2array(starfish_genes);

% only global coordinates (um) to pixel space-converted
% genes = nilsson_genes;
% N_withgenes=table(nilsson_x_sorted(iN_2), nilsson_y_sorted(iN_2), genes);
% genes = starfish_genes;
% S_withgenes=table(starfish_x_sorted_normalized(iS_2), starfish_y_sorted_normalized(iS_2), genes);

% only global coordinates (um) to pixel space-converted and rounded
% genes = nilsson_genes;
% N_withgenes=table(nilsson_x_sorted_rounded(iN_2), nilsson_y_sorted_rounded(iN_2), genes);
% genes = starfish_genes;
% S_withgenes=table(starfish_x_sorted_normalized_rounded(iS_2), starfish_y_sorted_normalized_rounded(iS_2), genes);


% % only global coordinates (um) to pixel space-converted and shifted
% genes = nilsson_genes;
% N_withgenes=table(nilsson_x_sorted(iN_2), nilsson_y_sorted(iN_2), genes);
% genes = starfish_genes;
% S_withgenes=table(starfish_x_sorted_normalized_shift(iS_2), starfish_y_sorted_normalized_shift(iS_2), genes);

% % 
% global coordinates (um) to pixel space-converted, shifted rounded
genes = nilsson_genes;
N_withgenes=table(nilsson_x_sorted_rounded(iN_2), nilsson_y_sorted_rounded(iN_2), genes);
genes = starfish_genes;
S_withgenes=table(starfish_x_sorted_normalized_shift_rounded(iS_2), starfish_y_sorted_normalized_shift_rounded(iS_2), genes);

genenames = unique(N_withgenes(:,3));

%nilsson as ground truth
S_withgenes_sub = S_withgenes(logical(ones(size(S_withgenes,1),1)),{'Var1', 'Var2'});

%starfish as ground truth
%N_withgenes_sub = N_withgenes(logical(ones(size(N_withgenes,1),1)),{'Var1', 'Var2'});

%% histogram plot function
percentageOverlap = cell(32,1); % used to save the percentage of perfectly overlapping reads
Overlaps = cell(32,1); 
perfectOverlaps = cell(32,1);
cmap=hsv(size(genenames,1));
figure()
for i = 1:size(genenames,1)
    index = strcmp(string(N_withgenes.genes),string(genenames.genes{i}));
    N_withgenes_sub = N_withgenes(index, {'Var1', 'Var2'});
    [C_withgenes_sub,iN_withgenes_sub,iS_withgenes_sub] = intersect(N_withgenes_sub, S_withgenes_sub,'rows');
    iS_withgenes_sub = find(ismember(S_withgenes_sub,C_withgenes_sub));
    counts = [];
    subplot(4,8,i)
    for j = 1:size(genenames,1)
        index = strcmp(string(S_withgenes.genes(iS_withgenes_sub)),string(genenames.genes{j}));
        counts = [counts; sum(index)];
    end
        
    subplot(4,8,i)
    hist(index, unique(index));
    bar([1:size(genenames)], counts,'b');
    hold on
    bar(i, counts(i),'r');
    title(string(genenames.genes{i}))
    xlabel('Starfish genes')
    ylabel('Overlapping reads')
    
    %sgtitle('Global coordinates (um) to pixel space-converted and rounded')
    %sgtitle('Global coordinates (um) to pixel space-converted and shifted')
    %sgtitle('Global coordinates (um) to pixel space-converted, shifted and rounded')

    sum(counts); 
    sum(counts(i));
    m(i,:) = [sum(counts)];
    m1(i,:) = [sum(counts(i))];
    percentagePerfectOverlap = [m1(i,:)/m(i,:)];
    percentagePerfectOverlap = percentagePerfectOverlap*100;
    percentagePerfectOverlap;
    overlaps{i} = m(i,:); 
    perfectOverlaps{i} = m1(i,:);
    overlaps = overlaps.';
    percentageOverlap{i} = percentagePerfectOverlap; 
    sumOverlapsAll = sum(m);
    sumPerfectOverlaps = sum(m1);
end 
    %% scatter plot function
cmap=hsv(size(genenames,1));
figure()
for i = 1:size(genenames,1)
    index = strcmp(string(N_withgenes.genes),string(genenames.genes{i}));
    N_withgenes_sub = N_withgenes(index, {'Var1', 'Var2'});
    [C_withgenes_sub,iN_withgenes_sub,iS_withgenes_sub] = intersect(N_withgenes_sub, S_withgenes_sub,'rows');
    iS_withgenes_sub = find(ismember(S_withgenes_sub,C_withgenes_sub));
    counts = [];
    subplot(4,8,i)
   
    gscatter(N_withgenes.Var1(index), N_withgenes.Var2(index), N_withgenes.genes(index),cmap(i,:), '.');
    title(string(genenames.genes{i}),'Color',cmap(i,:));
    hold on
    [tf, idx] = ismember(string(S_withgenes.genes(iS_withgenes_sub)),string(genenames.genes));
    gscatter(S_withgenes.Var1(iS_withgenes_sub), S_withgenes.Var2(iS_withgenes_sub), S_withgenes.genes(iS_withgenes_sub),cmap(idx,:),'o');
    
    indexxx = strcmp(string(S_withgenes.genes(iS_withgenes_sub)),string(genenames.genes{i}));
    gscatter(S_withgenes.Var1(iS_withgenes_sub(indexxx)), S_withgenes.Var2(iS_withgenes_sub(indexxx)), S_withgenes.genes(iS_withgenes_sub(indexxx)),cmap(i,:),'o');
    
      if i == 1
        x = xlim;
        y = ylim;
      end
      
    xlim([x(1) x(2)])
    ylim([y(1) y(2)]) 
    legend off;
    %sgtitle('Global coordinates (um) to pixel space-converted and rounded')
    %sgtitle('Global coordinates (um) to pixel space-converted and shifted')
    %sgtitle('Global coordinates (um) to pixel space-converted, shifted and rounded')

end

%% spot matching in numbers
% finding matching coordinates
[matchingCoordinates,iNN,iSS] = outerjoin(N_withgenes,S_withgenes,'LeftKeys',{'Var1','Var2'},'RightKeys',{'Var1','Var2'});

% finding spots that match in genes as well 
[matchingCoordinatesAndGenes,iNN_1,iSS_1] = innerjoin(N_withgenes,S_withgenes,'LeftKeys',{'Var1','Var2', 'genes'},'RightKeys',{'Var1','Var2', 'genes'});

size(matchingCoordinates)
size(matchingCoordinatesAndGenes)

%% Plotting both the histograms and the scatter plots
cmap=hsv(size(genenames,1));
figure()
for i = 1:size(genenames,1)
    
    %nilsson as ground truth
    index = strcmp(string(N_withgenes.genes),string(genenames.genes{i}));
    N_withgenes_sub = N_withgenes(index, {'Var1', 'Var2'});
    
    %starfish as ground truth
%     index = strcmp(string(S_withgenes.genes),string(genenames.genes{i}));
%     S_withgenes_sub = S_withgenes(index, {'Var1', 'Var2'});
        
    [C_withgenes_sub,iN_withgenes_sub,iS_withgenes_sub] = intersect(N_withgenes_sub, S_withgenes_sub,'rows');
    
    iS_withgenes_sub = find(ismember(S_withgenes_sub,C_withgenes_sub));
%     coords = [S_withgenes.Var1(iS_withgenes_sub) S_withgenes.Var2(iS_withgenes_sub)];
%     [u,a,b] = unique(coords, 'rows', 'first');
%     idxDupRows = setdiff(1:size(coords,1), a);
%     dupRowValues = unique(coords(idxDupRows,:), 'rows');
%     XXX = find(S_withgenes_sub.Var1 == dupRowValues(1) & S_withgenes_sub.Var2 == dupRowValues(2));
%     S_withgenes.genes(XXX);
%     
%     [uni, ~, idx] = unique(S_withgenes.genes(iS_withgenes_sub));
    counts = [];
    
    subplot(4,8,i)
    %plot(N_withgenes_sub.Var1, N_withgenes_sub.Var2, 'o', 'MarkerEdgeColor','black');
    
%     nilsson as ground truth
%     gscatter(N_withgenes_sub.Var1, N_withgenes_sub.Var2, N_withgenes.genes);
%     hold on
%     gscatter(S_withgenes.Var1(iS_withgenes_sub), S_withgenes.Var2(iS_withgenes_sub), S_withgenes.genes(iS_withgenes_sub));
       
%     gscatter(N_withgenes.Var1(index), N_withgenes.Var2(index), N_withgenes.genes(index),cmap(i,:), '.');
%     title(string(genenames.genes{i}),'Color',cmap(i,:));
%     hold on
%     [tf, idx] = ismember(string(S_withgenes.genes(iS_withgenes_sub)),string(genenames.genes));
%     gscatter(S_withgenes.Var1(iS_withgenes_sub), S_withgenes.Var2(iS_withgenes_sub), S_withgenes.genes(iS_withgenes_sub),cmap(idx,:),'o');
%     
% % gscatter plots alphabetical, i.e. co-ocurring genes are plotted as the last gene in the alphabet  
%     %-> simply plot the current gene on top of the duplicates
%     indexxx = strcmp(string(S_withgenes.genes(iS_withgenes_sub)),string(genenames.genes{i}));
%     gscatter(S_withgenes.Var1(iS_withgenes_sub(indexxx)), S_withgenes.Var2(iS_withgenes_sub(indexxx)), S_withgenes.genes(iS_withgenes_sub(indexxx)),cmap(i,:),'o');
%     
    %starfish as ground truth
%     gscatter(S_withgenes_sub.Var1, S_withgenes_sub.Var2);
%     hold on
%     gscatter(N_withgenes.Var1(iN_withgenes_sub), N_withgenes.Var2(iN_withgenes_sub), N_withgenes.genes(iN_withgenes_sub));

% commented out when using the barplot function and uncommented when using
% % %the scatter function
%     if i == 1
%         x = xlim;
%         y = ylim;
%     end
    
    for j = 1:size(genenames,1)
        index = strcmp(string(S_withgenes.genes(iS_withgenes_sub)),string(genenames.genes{j}));
        counts = [counts; sum(index)];
    end
    
    subplot(4,8,i)
    hist(index, unique(index));
    bar([1:size(genenames)], counts,'b');
    hold on
    bar(i, counts(i),'r');
    title(string(genenames.genes{i}))
    xlabel('Starfish genes')
    ylabel('Overlapping reads')
    sgtitle('Round to 1 digits to the left of the decimal point')
    
    sum(counts) 
    sum(counts(i))
    m(i,:) = [sum(counts)]
    m1(i,:) = [sum(counts(i))]
    sumOverlapsAll = sum(m)
    sumPerfectOverlaps = sum(m1)
    
% % commented out when using the barplot function and uncommented when using
% % the scatter function
%      xlim([x(1) x(2)])
%      ylim([y(1) y(2)]) 
%     legend off;
%     
end




