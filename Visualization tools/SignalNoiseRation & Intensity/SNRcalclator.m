% script for calculating the signal to noise ratio 
% intensity data generated using Xiaoyans scipt called 
% PrototypeSNRProfiler.m
close 
clear
for  s = 2:5 % four channels 
    for k = ["SBH","SBL"] % two sequencing chemistries 
        all_means = [];
        all_bg = [];
        all_SNR = [];
        tabledata = [];
        for r = ["ROI1", "ROI2", "ROI3"] % three regions of interest
                tabledata = readtable(strcat('/Volumes/SBHvsSBL/Orthogonal prpjection/', (k), '_Round6Cycle3_', (k), '_', (r), '/', (k), '_Round6Cycle3_', (k), '_', (r), '_c', num2str(s), '_ORG.csv')); 
                tabledata = cat(1, tabledata, tabledata);
%                 tabledata = rmoutliers(tabledata, 'grubbs');
                arraydata = table2array(tabledata, 'grubbs');
               
                % finding the mean signal
                signalmean = mean(arraydata, 1);
                
                % calculating the mean background
                backgroundPos1 = (signalmean (:,1));
                backgroundPos2 = (signalmean (:,2));
                backgroundPos3 = (signalmean (:,20));
                backgroundPos4 = (signalmean (:, 21));
                backgroundmean = ((backgroundPos1+backgroundPos2+backgroundPos3+backgroundPos4)/4);
                
                % calculating the signal to noise ratio 
                SNR = signalmean/backgroundmean;  
              
        end
     
        if k == "SBH" 
            SNR_SBH = SNR;
        else 
            SNR_SBL = SNR;
        end
        
    end  

    subplot(2,2, s-1);
    plot(SNR_SBH, 'MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red')
    hold on 
    plot(SNR_SBL, 'MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red')
    title(strcat('Signal to noise ratio ', ' channel  ', num2str(s)));
    legend SBH SBL
    xlabel ('Pixel number')
    ylabel ('Signal to noise ratio')
    ax = gca;
    ax.FontSize = 12; 
    xlim([1 21])
    ylim([1 9])
end 

%%