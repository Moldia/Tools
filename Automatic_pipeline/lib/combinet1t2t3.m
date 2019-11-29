t1=imread('L:\Processing hmHCA_10_sect2_big\Base_1_c1m1_ORG.tif');
t2=imread('L:\Processing hmHCA_10_sect2_big\Base_1_c1m2_ORG.tif');
t3=imread('L:\Processing hmHCA_10_sect2_big\Base_1_c1m3_ORG.tif');

overlap=0

xy1=size(t1)
xy2=size(t2)
xy3=size(t3)

x1=xy1(1);
x2=xy2(1);
x3=xy3(1);

startx2=floor(x1*(1-overlap));
startx3=floor(x2*(1-overlap));

shiftx1=x1-startx2

cutt1=t1(1:startx2,:);
cutt2=t2(1:startx3,:);

wholeimage=vertcat(cutt1,cutt2,t3);


%THIS PART IS WORK


%MODIFY IN DETAILS
genes1=readtable("L:\Processing hmHCA_10_sect2_big\decoding_1\decoding\QT_0.5_0_details.csv");
genes2=readtable("L:\Processing hmHCA_10_sect2_big\decoding_2\decoding\QT_0.5_0_details.csv");
genes3=readtable("L:\Processing hmHCA_10_sect2_big\decoding_3\decoding\QT_0.5_0_details.csv");
shiftx2=startx2;
shiftx3=startx2+startx3;
genes2.PosY=genes2.PosY+shiftx2;
genes3.PosY=genes3.PosY+shiftx3;
genesALL=vertcat(genes1,genes2,genes3);
writetable(genesALL,'L:\Processing hmHCA_10_sect2_big\decodingALL\QT_0.5_0_details.csv');


%PARAMETERS PRINT
outputdirectory='L:\Processing hmHCA_10_sect2_big\decodingALL'
q.quality_threshold = 0.5;        
q.general_stain_threshold = 0;
d.taglist = 'L:\Processing hmHCA_10_sect2_big\taglist_human120.csv';

 p.background_image = wholeimage;
        p.scale = 1;		% image scale
    p.I_want_to_plot_on_white_backgound = 1;
    % options
    p.exclude_NNNN_YN = 0;
    p.plot_reads_beforeQT_YN = 0;
    p.plot_ref_general_stain = 0; 

seqplotting2(outputdirectory, d.taglist, q, p)



