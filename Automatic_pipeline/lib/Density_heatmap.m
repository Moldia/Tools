[name,pos]=getinsitudata('E:\Whole_organoid_pseudoAnchor_v5\Decoding\QT_0.7_details_noNNNN.csv');
information=[name;pos(:,1);pos(:,2)];
N=10
xpos=pos(:,1);
ypos=pos(:,2);
x = linspace(0,1,N) ;
y = linspace(0,1,N) ;
T = rand(N) ;

name=linspace(1,1,size(xpos,1));
pcolor(xpos,ypos,name) ;

colorbar
pcolor(information(2,:),information(2,:),information(3,:)) ;
colorbar