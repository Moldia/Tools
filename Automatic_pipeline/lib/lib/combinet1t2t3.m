t1=imread('J:\MIPs_ms120\sergio\withanchor\Tile1_cyc_1_CH_1m1_ORG.tif');
t2=imread('J:\MIPs_ms120\sergio\withanchor\Tile2_cyc_1_CH_1m2_ORG.tif');
t3=imread('J:\MIPs_ms120\sergio\withanchor\Tile3_cyc_1_CH_1m3_ORG.tif');

overlap=0.05

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

%THIS IS NOT WORKING
imwrite(wholeimage,'J:\MIPs_ms120\sergio\CellprofilerALL\wholeimage.tif');
f=fopen('J:\MIPs_ms120\sergio\CellprofilerALL\wholeimage.tif','w');
fwrite(f,wholeimage,'uint16');
fclose(f);

%THIS PART IS WORK


%MODIFY IN DETAILS
genes1=readtable("J:\MIPs_ms120\sergio\Cellprofilerm1\QT_0.5_0_details.csv");
genes2=readtable("J:\MIPs_ms120\sergio\Cellprofilerm2\QT_0.5_0_details.csv");
genes3=readtable("J:\MIPs_ms120\sergio\Cellprofilerm3\QT_0.5_0_details.csv");
shiftx2=startx2;
shiftx3=startx2+startx3;
genes2.PosY=genes2.PosY+shiftx2;
genes3.PosY=genes3.PosY+shiftx3;
genesALL=vertcat(genes1,genes2,genes3);
writetable(genesALL,'J:\MIPs_ms120\sergio\CellprofilerALL\QT_0.5_0_details.csv');




%MODIFY NONA

genes1=readtable("J:\MIPs_ms120\sergio\Cellprofilerm1\QT_0.5_details_noNNNN.csv");
genes2=readtable("J:\MIPs_ms120\sergio\Cellprofilerm2\QT_0.5_details_noNNNN.csv");
genes3=readtable("J:\MIPs_ms120\sergio\Cellprofilerm3\QT_0.5_details_noNNNN.csv");
shiftx2=startx2;
shiftx3=startx2+startx3;
genes2.PosY=genes2.PosY+shiftx2;
genes3.PosY=genes3.PosY+shiftx3;
genesALL=vertcat(genes1,genes2,genes3);
writetable(genesALL,'J:\MIPs_ms120\sergio\CellprofilerALL\QT_0.5_details_noNNNN.csv');



%PARAMETERS PRINT
outputdirectory='J:\MIPs_ms120\sergio\CellprofilerALL'
q.quality_threshold = 0.5;        
q.general_stain_threshold = 0;
d.taglist = 'J:\MIPs_ms120\CellProfiler\taglist_mouse120.csv';

 p.background_image = wholeimage;
        p.scale = 1;		% image scale
    p.I_want_to_plot_on_white_backgound = 1;
    % options
    p.exclude_NNNN_YN = 0;
    p.plot_reads_beforeQT_YN = 0;
    p.plot_ref_general_stain = 0; 

seqplotting2(outputdirectory, d.taglist, q, p)



