% script for plotting the maximum intensities
% intensity data generated using Xiaoyans scipt called 
% PrototypeSNRProfiler.m
clear
close
for  s = 2:5 % four channels 
    for k = ["SBH","SBL"] % two sequencing chemistries 
        IntensityData = [];
        for r = ["ROI1", "ROI2", "ROI3"] % three regions of interest
                tabledata = readtable(strcat('/Volumes/SBHvsSBL/Orthogonal prpjection/', (k), '_Round6Cycle3_', (k), '_', (r), '/', (k), '_Round6Cycle3_', (k), '_', (r), '_c', num2str(s), '_ORG.csv')); 
                
                REMOVEDOUTLIERS = rmoutliers(tabledata, 'grubbs');

                middlePixelIntensity = REMOVEDOUTLIERS.Var11;
                %maximumIntensity = (middlePixelIntensity, 100);  
                
                IntensityData = cat(1, IntensityData, middlePixelIntensity);
                
                IntensityData = num2cell(IntensityData); 
                IntensityData = cell2mat(IntensityData);
        end
        
     
        if k == "SBH" 
            % preparing the data
            IntensityDataSBH = IntensityData;
            IntensityDataSBH_cell = num2cell(IntensityDataSBH);
            IntensityDataSBH_mat = cell2mat(IntensityDataSBH_cell);
            
            % preparing the cell strings
            length_IntensityData = length(IntensityData);
            cell_length_IntensityData = cell((length_IntensityData),1);
            names_SBH = strcat(cell_length_IntensityData, k, num2str(s));
            All_Cell_SBH = names_SBH;
        else 
            % preparing the data
            IntensityDataSBL = IntensityData; 
            IntensityDataSBL_cell = num2cell(IntensityDataSBL);
            IntensityDataSBL_mat = cell2mat(IntensityDataSBL_cell);
                       
            % preparing the cell strings
            length_IntensityData = length(IntensityData);
            cell_length_IntensityData = cell((length_IntensityData),1);
            names_SBL = strcat(cell_length_IntensityData, k, num2str(s));
            All_Cell_SBL = names_SBL;
        end
    
    end 
    
    All_Cell_SBHandSBL = [IntensityDataSBH_mat; IntensityDataSBL_mat]; 
    All_Cell_Names = [All_Cell_SBH; All_Cell_SBL];
    
    width = 0.1;
    showmean = true;
    violinalpha = 0.3;
    
    subplot(2,2, s-1);
    violinplot(All_Cell_SBHandSBL, All_Cell_Names);
    title(strcat('Maximum intensity of RCPs ', ' channel  ', num2str(s)));
    ylabel('Intensity');
    ylim([0, 5000]);

end
