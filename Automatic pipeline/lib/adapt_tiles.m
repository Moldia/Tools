function csvfile = adapt_tiles(input_directory,file_name,output_directory)
tiledata=readtable([input_directory,'\','tile_coordinates_',file_name,'.csv'])
ss=table((1:size(tiledata,1))',tiledata(:,1),tiledata(:,2))
ss.Properties.VariableNames = {'Metadata_position' 'Tile_xPos' 'Tile_yPos'};
ss=splitvars(ss);
ss.Properties.VariableNames = {'Metadata_position' 'Tile_xPos' 'Tile_yPos'};
csvfile=strcat(output_directory,'\tiled.csv')
writetable(ss,strcat(output_directory,'\tiled.csv'));
end
