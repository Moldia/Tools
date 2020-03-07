load E:\pciSeq\Data\AllSections\gSetCA1all.mat;
load F:\HCA_09b_mouse_120genes\pciseq\Basecalling_NewRegistration2\oThere.mat
load('F:\HCA_09b_mouse_120genes\pciseq\Basecalling_NewRegistration2\Cellcalling\CellMap_left.mat')


gSet.GeneExp=readmatrix('G:\Glioblastoma_data\Sequencing data/clusters_expression_Sample1.csv');
gSet.GeneExp;
gSet.GeneExp=gSet.GeneExp';
gSet.GeneExp=gSet.GeneExp(2:end,:);
gSet.GeneExp=gSet.GeneExp';
gSet.nGenes=size(gSet.GeneExp,1);
gSet.nCells=size(gSet.GeneExp,2);
%extrainfo=readtable('C:\Users\sergio.salas\Downloads\sample-annotations\sample_annotations.csv');
%gSet.CellInfo=[]
GENAME=readtable('G:\Glioblastoma_data\Sequencing data/clusters_expression_Sample1.csv');
gSet.GeneName=[GENAME.Var1];
gSet.Class=GENAME.Properties.VariableNames;
%CENAME=readcell('F:/Mouse_Brain_cluster_sc1.csv');
 %LEVEL 1 
gSet.Class=gSet.Class(2:end);
ClassNames=unique(gSet.Class)';
gSet.CellName=gSet.CellName(1:size(gSet.Class,2))';
gKeep=gSet;


%LEVEL1
minimrat=1;
minimexp=1;
maxexp=0.05;
gSetSt=gSet;
%gSet.Class=CENAME(2:end,2)'; %LEVEL 1 needs to be 2
ClassNames=unique(gSet.Class)';
CLSEMPTY=zeros(1,size(gSet.Class,2));
[CLS1,FINAL,TYPE,gSet,gStable] = gene_selection2(gSet,ClassNames,o,gKeep,CLSEMPTY,minimrat,minimexp,maxexp);
path='GLIOBLASTOMA_SAMPLE1';
name='GLIOBLASTOMA_SAMPLE1';
[FINAL1]=gene_curation(gSet,gStable,FINAL1,path,name,ClassNames,o,CellMap);


%LEVEL2
gSet.Class=CENAME(2:end,3)'; %LEVEL 2 
ClassNames=unique(gSet.Class)';
%CLS1=CENAME(2,2:end); %ALWAYS ASSIGN ONE COLUMN LESS
[CLS2,FINAL,TYPE,gSet,gStable] = gene_selection2(gSet,ClassNames,o,gKeep,CLS1,minimrat,minimexp,maxexp);
path='LEV1_MOUSEB';
name='';
[FINAL2]=gene_curation(gSet,gStable,FINAL,path,name,ClassNames,o,CellMap);



%LEVEL3
gSet.Class=CENAME(2:end,4)'; %LEVEL 3
ClassNames=unique(gSet.Class)';
%CLS1=CENAME(2,2:end); %ALWAYS ASSIGN ONE COLUMN LESS
[CLS3,FINAL,TYPE,gSet,gStable] = gene_selection2(gSet,ClassNames,o,gKeep,CLS2,minimrat,minimexp,maxexp);
path='LEV1_MOUSEB';
name='';
[FINAL3]=gene_curation(gSet,gStable,FINAL,path,name,ClassNames,o,CellMap);



%LEVEL4
gSet.Class=CENAME(2:end,5)'; %LEVEL 4
ClassNames=unique(gSet.Class)';
%CLS1=CENAME(2,2:end); %ALWAYS ASSIGN ONE COLUMN LESS
[CLS4,FINAL,TYPE,gSet,gStable] =  gene_selection2(gSet,ClassNames,o,gKeep,CLS3,minimrat,minimexp,maxexp);
path='LEV1_MOUSEB';
name='';
[FINAL4]=gene_curation(gSet,gStable,FINAL,path,name,ClassNames,o,CellMap);



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

