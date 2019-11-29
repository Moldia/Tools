comasep = readtable('I:\3x3tile\Preprocess\withANCHOR\positions4.csv');
comasep;
im = imread('I:\3x3tile\190919_msSBH_02_rd7cy2-Orthogonal Projection-01\190919_msSBH_02_rd7cy2-Orthogonal Projection-01_c1+2+3+4+5.tif');

%imr=imresize(im,7);
image(im);

hold on
comasep.gene=categorical(comasep.gene);
gscatter(comasep.truex,comasep.truey,comasep.gene);
hold off

class(comasep.gene)