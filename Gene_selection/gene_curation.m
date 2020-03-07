
function[FINAL,SIM,EXP]=gene_curation(gSet,gStable,FINAL,path,name,ClassNames,o,CellMap)
ExcludeGenes={};

%CL=ClassNames(60:88);

disp(strcat("Now we use the Desired genes only"));
gSet=gStable;

EXGENE=(~ismember(gSet.GeneName,ExcludeGenes));
gSet.GeneName=gSet.GeneName(EXGENE,:);
gSet.GeneExp=gSet.GeneExp(EXGENE,:);
gSet.nGenes=size(gSet.GeneExp,1);
gSet.nCells=size(gSet.GeneExp,2);
GeneNames=unique(gSet.GeneName);
nG = length(GeneNames);
nK = length(ClassNames); % last is zero-expression
ClassNames=unique(gSet.Class)';
CL=ClassNames(:);
%%%%%%%%%%%%THIS WOULD BE THE POINT TO START THE LOOP%%%%%%%%%%%%%%%%%%%%

GeneNames=[];
for col=1:size(FINAL,2)
   GeneNames=[GeneNames,FINAL(:,col)']; 
end
ClassPrior = [.5*ones(1,nK-1)/nK .5];
ClassDisplayNames = ClassNames;
%In case of adding clustering data, add here
gSet.GeneExp=gSet.GeneExp;
gSub = gSet.GeneSubset(GeneNames);
gSub.GeneExp=gSub.GeneExp';
g=gSub;
%THIS IS MEAN CALCULATION
p=1;
q=0;
nG= size(gSub.GeneName,1);
nK = length(ClassNames); % last is zero-expression
%nC = size(CellYX,1)+1; % last is misreads
%nS = size(SpotYX,1);
MeanClassExp = zeros(nK, nG);
o.Inefficiency=0.2;
for k=1:nK-1 % don't include last since it is zero-expression class
            Cells = g.IdentifyCells(ClassNames{k});
            h = g;
            h.GeneExp = g.GeneExp(Cells,:);
            h.GeneName = g.GeneName;
            h.nGenes = g.nGenes;
            h.nCells = length(Cells);
            if ~isempty(g.CellName), h.CellName = g.CellName(Cells); end
            if ~isempty(g.tSNE), h.tSNE = g.tSNE(:,Cells); end

            h.GeneExp=h.GeneExp';
            Norm = mean(h.GeneExp.^q,1).^(1/q);
            Not0Norm = (Norm>0);
            s = h.CellSubset(Not0Norm);
           % s.GeneExp=s.GeneExp';
            s.GeneExp=s.GeneExp';
            s.GeneExp = bsxfun(@rdivide, h.GeneExp(:,Not0Norm), Norm(Not0Norm).^p);
            ScaleFac = (sum(h.GeneExp(:).^q)/sum(s.GeneExp(:).^q)).^(1/q);
            s.GeneExp = s.GeneExp*ScaleFac;
    MeanClassExp(k,:) = o.Inefficiency * mean(s.GeneExp,2)';
    
end

lMeanClassExp = log(MeanClassExp + o.SpotReg); 


% Generation of random cells spread around the sample
simulatedXY=rand(1000,2);
simulatedXY=simulatedXY*1000;
simulatedCT=[];
for s=1:1000
%selection=randsample(ClassNames,1);
selection=randsample(CL,1); %This should change depending on the cycle
simulatedCT=[simulatedCT;selection];
end
disp(strcat('Selection was ',selection));
%% Initiate variables for the loop

TOTgene=[];
TOTpos=[];
%Now it's time to calculate the expression on each point
for j=1:size(MeanClassExp,2)
disp(j)
for x =1:size(MeanClassExp,1)
class=ClassNames(x);
v=find(strcmp(simulatedCT, class));
% disp(j)
for num=1:size(v)
   
    cellsel=v(num);
    lambda=MeanClassExp(x,j);
    r = poissrnd(lambda);
    if r>0
    for selection= 1:r
        TOTgene=[TOTgene;GeneNames(j)];
        TOTpos=[TOTpos;simulatedXY(cellsel,1)+((rand(1)-0.5)*10), simulatedXY(cellsel,2)+((rand(1)-0.5)*10)];
    end
    end 
end
end
end

%% This plot represents the simulated dataset

% figure(12498);
% scatter(TOTpos(:,2),TOTpos(:,1),'b');
% hold on
% scatter(simulatedXY(:,2),simulatedXY(:,1),'r');

%Simulated dataset is saved in the appropiate variables for performing
%pciSeq
CellYX=simulatedXY;
SpotYX=TOTpos;
SpotGeneName=TOTgene;


%% NOW PERFORMING PCISEQ


rp = regionprops(CellMap);
CellArea0 = vertcat(rp.Area);
CellArea0=CellArea0(1:size(CellYX,1));

MeanCellRadius = mean(sqrt(CellArea0/pi))*.5; % the dapi part is only half of the typical radius
RelCellRadius = [sqrt(CellArea0/pi)/MeanCellRadius; 1]; % but here we want the whole thing

%% get arrays ready

% SpotGene(nS): which gene is each spot
% MeanClassExp(nK,nG): mean expression of each gene in each class
% Neighbors(nS, nN): closest neighboring cells for each spot
% D(nS, nN): distance penalty for each of these
% GeneNames(nG): name of each gene
% ClassNames(nK): name of each class

[GeneNames, ~, SpotGeneNo] = unique(SpotGeneName);
TotGeneSpots = accumarray(SpotGeneNo,1);
ClassNames = vertcat(unique(ClassNames, 'stable'), {'Zero'});

gSub = gSet.GeneSubset(GeneNames);
g=gSub;
%TotGeneSpots=[TotGeneSpots;1] % I added this ine
nG = length(unique(gSub.GeneName));
nK = length(ClassNames); % last is zero-expression
nC = size(CellYX,1)+1; % last is misreads
nS = size(SpotYX,1);
nN = o.nNeighbors+1; % last is misreads (always a neighbor)
ClassPrior = [.5*ones(1,nK-1)/nK .5];
ClassDisplayNames = ClassNames;
MeanClassExp = zeros(nK, nG);

MeanClassExp = zeros(nK, nG);
for k=1:nK-1 % don't include last since it is zero-expression class
            Cells = g.IdentifyCells(ClassNames{k});
            h = g;
            h.GeneExp = g.GeneExp(:,Cells);
            h.GeneName = g.GeneName;
            h.nGenes = g.nGenes;
            h.nCells = length(Cells);
            if ~isempty(g.CellName), h.CellName = g.CellName(Cells); end
            if ~isempty(g.tSNE), h.tSNE = g.tSNE(:,Cells); end

           % h.GeneExp=h.GeneExp';
            Norm = mean(h.GeneExp.^q,1).^(1/q);
            Not0Norm = (Norm>0);
            s = h.CellSubset(Not0Norm);
           % s.GeneExp=s.GeneExp';
            s.GeneExp=s.GeneExp';
            s.GeneExp = bsxfun(@rdivide, h.GeneExp(:,Not0Norm), Norm(Not0Norm).^p);
            ScaleFac = (sum(h.GeneExp(:).^q)/sum(s.GeneExp(:).^q)).^(1/q);
            s.GeneExp = s.GeneExp*ScaleFac;
            TF = isempty(s.GeneExp);
            if TF==1;
                s.GeneExp=zeros(size(s.GeneExp,1),1);
            end
    MeanClassExp(k,:) = o.Inefficiency * mean(s.GeneExp,2)';
    
end










lMeanClassExp = log(MeanClassExp + o.SpotReg); 

% now find each spot's neighboring cells and distances (nS, nN)
[Neighbors, Dist] = knnsearch(CellYX, SpotYX, 'K', nN);
Neighbors(:,end) = nC; % set last neighbor to misreads

D = -Dist.^2./(2*MeanCellRadius^2) - log(2*pi*MeanCellRadius^2); % don't normalize: bigger cells express more
D(:,end) = log(o.MisreadDensity); % this is log likelihood of misread

% any inside cell radius given a bonus
%SpotInCell = IndexArrayNan(CellMap, (SpotYX - [round(y0) round(x0)])');
SpotInCell=zeros(size(SpotYX,1),1);
if Neighbors(SpotInCell>0,1)~=SpotInCell(SpotInCell>0)
    error('a spot is in a cell not closest neighbor!');
end
D(SpotInCell>0, 1) = D(SpotInCell>0, 1) + o.InsideCellBonus;

LogClassPrior = log(ClassPrior);

% this is area factor relative to that of the average cell
CellAreaFactor = (exp(-RelCellRadius.^2/2)*(1-exp(o.InsideCellBonus)) + exp(o.InsideCellBonus)) ...
    / (exp(-1/2)*(1-exp(o.InsideCellBonus)) + exp(o.InsideCellBonus));



%% initialize variables for main loop

pSpotNeighb = zeros(nS, nN); % prob each spot goes to each neighboring cell: last assigned to noise
pCellClass = zeros(nC, nK); % prob each cell goes to each class: last has zero expression

% start a spot in cell it is in, otherwise misread
pSpotNeighb(Neighbors==SpotInCell)=1;
pSpotNeighb(SpotInCell==0,end)=1;

% gammas start off as priors
eSpotGamma = ones(nC, nK, nG);
elSpotGamma = ones(nC, nK, nG)*psi(1); % start with r=1 prior, no spots

eGeneGamma = ones(nG,1); % start with just 1

% this is to check convergence
pSpotNeighbOld = zeros(nS, nN);

%% now main loop
for i=1:o.CellCallMaxIter
    % CellGeneCount(nC, nG): number of copies of each gene in each cell
    CellGeneCount = zeros(nC,nG);
    for n=1:nN-1
        c = Neighbors(:,n);
        CellGeneCount = CellGeneCount + accumarray([c, SpotGeneNo], pSpotNeighb(:,n), [nC,nG]);
    end
    
    %% call cell gammas

    % eSpotGamma(nC, nK, nG); expected gamma parameter
    % elSpotGamma(nC, nK, nG); expected log gamma parameter
    ScaledMean = CellAreaFactor.*reshape(MeanClassExp,[1 nK nG]);
    eSpotGamma = (o.rSpot+reshape(CellGeneCount,[nC 1 nG]))./(o.rSpot + ScaledMean);
    elSpotGamma = psi(o.rSpot+reshape(CellGeneCount,[nC 1 nG])) - log(o.rSpot + ScaledMean); % expectation of log gamma

    %% call cells
    
    % ScaledExp(nC, nK, nG): expected expression under current parameters
    ScaledExp = reshape(MeanClassExp,[1 nK nG]) .* reshape(eGeneGamma,[1 1 nG]) .* CellAreaFactor   + o.SpotReg;

    % pNegBin(nC, nK, nG): negbin parameter
    pNegBin = ScaledExp ./ (o.rSpot + ScaledExp);
    
    % wCellClass(nC, nK): summed log likelihoods
    wCellClass = sum(reshape(CellGeneCount,[nC 1 nG]).*log(pNegBin) + o.rSpot*log(1-pNegBin),3) + LogClassPrior;
    
    % pCellClass(nC, nK): probabilities
    pCellClass = LogLtoP(wCellClass')';

    
    
    %% call spots
    % wSpotCell(nS, nN)
    aSpotCell = zeros(nS, nN);
    for n=1:nN-1 % don't include misread possibility
        c = Neighbors(:,n);
        aSpotCell(:,n) = sum(pCellClass(c,:) .* lMeanClassExp(:,SpotGeneNo)',2) + ...
            sum(pCellClass(c,:) .* bi(elSpotGamma, c, 1:nK, SpotGeneNo), 2);
    end
    wSpotCell = aSpotCell + D ;
    
    pSpotNeighb = LogLtoP(wSpotCell')';
    MeanProbChanged = max(abs(pSpotNeighb(:)-pSpotNeighbOld(:)));
    fprintf('Iteration %d, mean prob change %f\n', i, MeanProbChanged)
    Converged = ( MeanProbChanged<o.CellCallTolerance);
    pSpotNeighbOld = pSpotNeighb;
    
        %% call gene gammas (etas)
    % to count non-background expression of each gene first compute background
    TotPredictedB = accumarray(SpotGeneNo, pSpotNeighb(:,end), [nG 1]);
    % and total spots in zero cells:
    pCellZero = pCellClass(:,nK); % prob a cell is class zero (nC)
    pSpotZero = sum(pSpotNeighb(:,1:nN-1).*pCellZero(Neighbors(:,1:nN-1)),2); % prob a spot comes from cell class zero (nS)
    TotPredictedZ = accumarray(SpotGeneNo, pSpotZero);
    
    % total counts predicted by all cells of each class (nK, nG)
    ClassTotPredicted = shiftdim(sum(eSpotGamma.*pCellClass.*CellAreaFactor,1),1).*(MeanClassExp + o.SpotReg);
    % total of each gene (nG): 
    TotPredicted = sum(ClassTotPredicted(1:nK-1,:),1)';
    
    eGeneGamma = (o.rGene + TotGeneSpots - TotPredictedB )./(o.rGene + TotPredicted); % First term was -TotPredictedZ
    if 0
        for gg=1:nG
            fprintf('%s:\t%f\n', GeneNames{gg}, eGeneGamma(gg)); 
        end
    end
    
    %% diagnostics
    if ~isempty(o.CellCallShowCenter) && (Converged || o.Graphics==2 || i==o.CellCallMaxIter)
        figure(3985471)
        
        % lines to show spots connected to cells
        [~, BestNeighb] = max(pSpotNeighb,[],2);
        SpotBestNeighb = bi(Neighbors,(1:nS)',BestNeighb(:));
        rn = SpotBestNeighb<nC & max(abs(SpotYX-o.CellCallShowCenter),[],2)<o.CellCallShowRad;
    end
    
    if Converged; break; end
end



o.pSpotCell = sparse(repmat(1:nS,1,nN)', Neighbors(:), pSpotNeighb(:));
o.CellYX = CellYX;
o.pCellClass = pCellClass;
o.ClassNames = ClassNames;




%%Here we start the analysis
MyCell=50;
Real_predicted= (o.rGene + TotGeneSpots)./(o.rGene + TotPredicted);
celltypes=zeros(size(o.pCellClass,1),1);
celltypes2=zeros(size(o.pCellClass,1),size(pNegBin,2)-1);
WeightALL=zeros(size(o.pCellClass,1),size(pNegBin,3));
WeightALL2=zeros(size(o.pCellClass,1),size(pNegBin,3),size(pNegBin,2)-1);


for dap = 1:size(o.pCellClass,1)
   celltype=find(o.pCellClass(dap,:)==max(o.pCellClass(dap,:)));
   celltype2=find(o.pCellClass(dap,:)~=max(o.pCellClass(dap,:))); 
   Myp = squeeze(pNegBin(dap,:,:));
   WeightMap = CellGeneCount(dap,:) .* log(Myp) +  o.rSpot*log(1-Myp);
   Weights=WeightMap(celltype,:)';
   Weights2=WeightMap(celltype2,:)';
   celltypes(dap,1)=celltype;
   celltypes2(dap,:)=celltype2;
   WeightALL(dap,:)=Weights;
   WeightALL2(dap,:,:)=Weights2;
end

noZeros=find(celltypes<max(celltypes));
    
Rcelltypes=celltypes(noZeros,:);
Rcelltypes2=celltypes2(noZeros,:);
RWeightALL=WeightALL(noZeros,:);
RWeightALL2=WeightALL2(noZeros,:,:);

RWM=mean(RWeightALL,1);
RWM2sub=mean(RWeightALL2,3);
RWM2=mean(RWM2sub,1);

 
 Diff=(abs(RWM)-abs(RWM2));   
 
 for i=1:size(RWM,1)
     for j=1:size(RWM,2)
    Diff(i,j)=Diff(i,j)/RWM(i,j);     
    end
 end

 
 celly=unique(Rcelltypes);
 CGWM=zeros(max(Rcelltypes),size(RWeightALL,2));
 for var=1:size(celly,1);
      vari=celly(var);
      sub=RWeightALL(find(Rcelltypes==vari),:);
      sub2=RWM2sub(find(Rcelltypes==vari),:);
      subM=mean(sub,1);
      subM2=mean(sub2,1);
      subTM=(subM-subM2);
      for j=1:size(subTM)
         subTM(1,j)=subTM(1,j)/subM(1,j); 
          
      end
      CGWM(vari,:)=subTM;
 end 

 
 % Selecting genes maximizing differences between class and others
 figure(121);
            imagesc(CGWM(celly,:));
            set(gca, 'xtick', 1:nG); set(gca, 'XTickLabel', GeneNames); set(gca, 'XTickLabelRotation', 90);
            set(gca, 'ytick', 1:nK); set(gca, 'yTickLabel', ClassNames(celly,:));
            title(sprintf('Cell %d: Contribution of genes to class scores', MyCell));


 
 
 
 
 
 
 
 
 
 
 
%top 10 genes
mostneeded=zeros(size(CGWM,1),size(CGWM,2));

for k = 1:size(celly,1);
vari=celly(k);
pos=[];
ve=sort(CGWM(vari,:));
fi=ve(end-5:end);
pos=[pos,find(CGWM(vari,:)==fi(end))]; 
pos=[pos,find(CGWM(vari,:)==fi(end-1))];
pos=[pos,find(CGWM(vari,:)==fi(end-2))];
pos=[pos,find(CGWM(vari,:)==fi(end-3))];
%pos=[pos,find(CGWM(vari,:)==fi(end-4))];
%pos=[pos,find(CGWM(vari,:)==fi(end-5))];


mostneeded(vari,pos(1))=20;
 mostneeded(vari,pos(2))=19;
 mostneeded(vari,pos(3))=18;
 mostneeded(vari,pos(4))=17;
% mostneeded(vari,pos(5))=16; 
%mostneeded(vari,pos(6))=15;

end
            
figure(12331);
            imagesc(mostneeded(celly,:));
            set(gca, 'xtick', 1:nG); set(gca, 'XTickLabel', GeneNames); set(gca, 'XTickLabelRotation', 90);
            set(gca, 'ytick', 1:nK); set(gca, 'yTickLabel', ClassNames(celly,:));
            title(sprintf('Cell %d: Most important genes on each class', MyCell));

% FINALGENES=GeneNames(importants);
% FINAL=[FINAL,FINALGENES];
% AMOUNT=sum(celltypes==classnum);
% CNT=[CNT,AMOUNT];



se=sum(mostneeded,1);
exclude=se<1;
ExcludeGenes=GeneNames(exclude);

assign={};
assign.CellYX=o.CellYX;
assign.pCellClass = o.pCellClass;
assign.ClassNames = o.ClassNames;
csvwrite(['F:\HCA_09b_mouse_120genes\pCellClass_',name,'.csv'],o.pCellClass);
T=cell2table(o.ClassNames);
writetable(T,['F:\HCA_09b_mouse_120genes\ClassNames_',name,'.csv']);
csvwrite(['F:\HCA_09b_mouse_120genes\CellYX_SELECT_',name,'.csv'],o.CellYX);
%save(assign.CellYX, assign.pCellClass,assign.ClassNames);



writecell(simulatedCT,['F:\HCA_09b_mouse_120genes\SimulatedCT_',name,'.csv']);

PREDICTED=[];
for w=1:size(o.pCellClass,1)
    PREDICTED=[PREDICTED,o.ClassNames(find(o.pCellClass(w,:)==max(o.pCellClass(w,:))),:)];
    
end


PREDICTED=PREDICTED(:,1:(end-1));

unipred=unique(PREDICTED);
unisimu=unique(simulatedCT);
results=zeros(size(unique(simulatedCT),1), size(unique(simulatedCT),1));

for i=1:size(unisimu,1)
  simu=unisimu(i);
  simumem=ismember(simulatedCT,simu);
  predic=PREDICTED(:,simumem);
  for j=1:size(unisimu,1)  
    results(i,j)=sum(ismember(predic,unisimu(j))); 
  end
end

SIM=[];
EXP=[];
for i=1:size(unisimu,1)
    for j=1:size(unisimu,1)
        if i~=j & results(i,j)~=0
            SIM=[SIM,unisimu(i)];
            EXP=[EXP,unisimu(j)];
        end
    end
end


FINAL=FINAL(~ismember(FINAL,ExcludeGenes));


end