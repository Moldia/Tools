function exclude_overlapping_regions(pathtotiled,tiledcsv,pathtoblobs,blobscsv)
tilespos=readtable([pathtotiled, tiledcsv]);
blobs=readtable([pathtoblobs,blobscsv]);
tilespos.Tile_xPos=round(tilespos.Tile_xPos,-2);
tilespos.Tile_yPos=round(tilespos.Tile_yPos,-2);

%for elements in the down corners, select them
el=unique(tilespos.Tile_xPos)
ydown=[]
for xposition=1:size(el)
   els=el(xposition)
   subsetx=tilespos(tilespos.Tile_xPos==els,:);
   edgex=subsetx(subsetx.Tile_yPos==max(subsetx.Tile_yPos),:); 
   ydown=[ydown,edgex.Metadata_position]
end

%For elements in the right top, select them
ele=unique(tilespos.Tile_yPos)
xright=[]
for yposition=1:size(ele)
   els=ele(yposition)
   subsety=tilespos(tilespos.Tile_yPos==els,:);
   edgey=subsety(subsety.Tile_xPos==max(subsety.Tile_xPos),:); 
   xright=[xright,edgey.Metadata_position];
end


tilenumber=unique(blobs.Metadata_position);
% Just select those points that are not in more than one tile, by removing
% the bottom and right corner of images that lay in the middle
blobs2=blobs(blobs.Location_Center_X<1894 | ismember(blobs.Metadata_position,xright),:);
blobs3=blobs2(blobs2.Location_Center_Y<1894 | ismember(blobs2.Metadata_position,ydown),:);
writetable(blobs3,[pathtoblobs,'blobsNO.csv'])

end