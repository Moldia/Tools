%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%this script allows to select genes%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% from single cell data and their clusters%%%%%%

%CREATE EMPTY OBJECTS. We need any random Cellmap for getting cell sizes
load F:\HCA_09b_mouse_120genes\pciseq\Basecalling_NewRegistration2\Cellmap.mat
o=iss;
gSet=GeneSet;

gSet.GeneExp=readmatrix('F:/Mouse_Brain_cluster_sc1.csv');
gSet.GeneExp;
gSet.GeneExp=gSet.GeneExp';
gSet.GeneExp=gSet.GeneExp(2:end,:);
gSet.GeneExp=gSet.GeneExp';
gSet.nGenes=size(gSet.GeneExp,1);
gSet.nCells=size(gSet.GeneExp,2);
%extrainfo=readtable('C:\Users\sergio.salas\Downloads\sample-annotations\sample_annotations.csv');
%gSet.CellInfo=[]
GENAME=readtable('F:/Mouse_Brain_genames1.csv');
gSet.GeneName=[GENAME.V1];
%CENAME=readcell('F:/Mouse_Brain_cluster_sc1.csv');
CENAME=readcell('F:/Taxonomy_rank_Mouse_brain.csv');
gSet.Class=CENAME(2:end,2)'; %LEVEL 1 
gSet.CellName=gSet.Class(1:size(gSet.Class,2))';
ClassNames=unique(gSet.Class)';
gKeep=gSet;

%%%%Basically we have gene selection for each depth in the cluster tree.
% We put on gSet.Class the class for each cell on the desired level
%FINAL1 include the genes curated for this level. We repreat this structure
% for each level we have

%LEVEL1
gSetSt=gSet;
gSet.Class=CENAME(2:end,6)'; %LEVEL 1 needs to be 2
ClassNames=unique(gSet.Class)';
CLSEMPTY=zeros(1,size(gSet.Class,2));
[CLS1,FINAL,TYPE,gSet,gStable] = gene_selectionHB(gSet,ClassNames,o,CellMap,gKeep,CLSEMPTY);
path='LEV1_MOUSEB';
name='';
[FINAL1]=gene_curation(gSet,gStable,FINAL,path,name,ClassNames,o,CellMap);


%LEVEL2
gSet.Class=CENAME(2:end,3)'; %LEVEL 2 
ClassNames=unique(gSet.Class)';
%CLS1=CENAME(2,2:end); %ALWAYS ASSIGN ONE COLUMN LESS
[CLS2,FINAL,TYPE,gSet,gStable] = gene_selection(gSet,ClassNames,o,CellMap,gKeep,CLS1);
path='LEV1_MOUSEB';
name='';
[FINAL2]=gene_curation(gSet,gStable,FINAL,path,name,ClassNames,o,CellMap);



%LEVEL3
gSet.Class=CENAME(2:end,4)'; %LEVEL 3
ClassNames=unique(gSet.Class)';
%CLS1=CENAME(2,2:end); %ALWAYS ASSIGN ONE COLUMN LESS
[CLS3,FINAL,TYPE,gSet,gStable] = gene_selection(gSet,ClassNames,o,CellMap,gKeep,CLS2);
path='LEV1_MOUSEB';
name='';
[FINAL3]=gene_curation(gSet,gStable,FINAL,path,name,ClassNames,o,CellMap);



%LEVEL4
gSet.Class=CENAME(2:end,5)'; %LEVEL 3
ClassNames=unique(gSet.Class)';
%CLS1=CENAME(2,2:end); %ALWAYS ASSIGN ONE COLUMN LESS
[CLS4,FINAL,TYPE,gSet,gStable] = gene_selection(gSet,ClassNames,o,CellMap,gKeep,CLS3);
path='LEV1_MOUSEB';
name='';
[FINAL4]=gene_curation(gSet,gStable,FINAL,path,name,ClassNames,o,CellMap);



%LEVEL5
gSet.Class=CENAME(2:end,6)'; %LEVEL 3
ClassNames=unique(gSet.Class)';
%CLS1=CENAME(2,2:end); %ALWAYS ASSIGN ONE COLUMN LESS
[CLS5,FINAL,TYPE,gSet,gStable] = gene_selection(gSet,ClassNames,o,CellMap,gKeep,CLS4);
path='LEV1_MOUSEB';
name='';
[FINAL5]=gene_curation(gSet,gStable,FINAL,path,name,ClassNames,o,CellMap);


%%%%%%%%%%%%%%%%%%%%ONCE WE HAVE IT FOR ALL THE LEVELS, WE MERGE ALL GENES
%%%%%%%%%%%%%%%%%%%%SELECTED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%AND PREDICT ON THE LAST LEVEL  %%%%%%%%%%%%%%%%%%%%%%%%


%FINAL LEVEL
FINALT=[FINAL1,FINAL2,FINAL3,FINAL4];
%FINALT=FINALT(~ismember(FINALT,exclude));
path='FINAL_LUNG';
name='FINAL_LUNG';

[FINALT2]=gene_curation(gSet,gStable,FINALT,path,name,ClassNames,o,CellMap);
mem=ismember(FINALT,FINALV1);
FINALCOM=FINALT(mem);
sum(ismember(unique(FINAL1),FINALCOM))
sum(ismember(unique(FINAL2),FINALCOM))
sum(ismember(unique(FINAL3),FINALCOM))

sum(ismember(FINALCOM,unique(FINAL3))& ismember(FINALCOM,unique(FINAL2)))
sum(ismember(FINALCOM,unique(FINAL1))& ismember(FINALCOM,unique(FINAL2)))
sum(ismember(FINALCOM,unique(FINAL1))& ismember(FINALCOM,unique(FINAL3)))
sum(ismember(FINALCOM,unique(FINAL1))& ismember(FINALCOM,unique(FINAL3))& ismember(FINALCOM,unique(FINAL2)))


[FINALT2]=gene_curation(gSet,gStable,FINALALL,path,name,ClassNames,o,CellMap);

gSub=gSet.GeneSubset(FINALT);
max(max(gSub.GeneExp'*0.02))

clustergram(MeanClassExp,'RowLabels',ClassNames,'ColumnLabels',GeneNames)
colormap(parula)






