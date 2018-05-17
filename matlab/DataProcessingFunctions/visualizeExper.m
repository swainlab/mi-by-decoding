function [ allCellsHeatMap,indexes ] = visualizeExper( infdata, xLims,yLims,nSingles ,TFname,classesLegend)
%VISUALIZEEXPER Visualization of single cell trajectories. 
% The function takes the infdata structure and plots a sample of single
% cell for each experiment + heatmap with all cells for each experiment. 
% nSingles: how many single cells to plot
%
% Author: Alejandro Granados
% Copyright Gasper Tkacik and Peter Swain 2018

nRows = 2;  
tp =2.5; %2.5 minutes timepoint (for plotting in minutes)
pCol = cool; %colormap
cc = pCol(ceil(rand(nSingles,1)*64),:); %sample colors
allCellsHeatMap =[]; 

idxAll=1; 
figure(1)
strains= fieldnames(infdata); %how many experiments in the structure

    for i = 1:length(strains)

        %concatenate the two conditions in a single matrix for plotting
        varSync = [infdata.(strains{i}).rich,infdata.(strains{i}).stress]; 
        %time axis: 
        midtimes = -(size(infdata.rep1.rich,2)-1): size(infdata.rep1.stress,2); 

        varMean = nanmean(varSync,1);
        
        %optional variables for shaded area plot: 
        % varN = sum(~isnan(varSync),1);       
        % varSErr = nanstd(varSync,0,1)./sqrt(varN);
        % varStd = nanstd(varSync,0,1); 
        
        subplot(nRows,length(strains),i);

        title(classesLegend{i})
        %optional: plot standard deviation of (stdErr) as shaded area
        %handles(i) = boundedline(midtimes,varMean,varStd,'black','alpha');
        hold on ; 
     
        %sample random single cell from dataset and plot them
            for pp = 1:nSingles
               idx=rand(1)*size(varSync,1); 
               thisCell = varSync(ceil(idx),:); 
               plot(midtimes,thisCell,'Color', [cc(pp,:),0.8],'LineWidth',0.5); 
            end
            
        plot(midtimes,varMean,'Color','black','LineWidth',1.5) ; 
        xlim([xLims(1),xLims(2)])
        ylim([yLims(1),yLims(2)])
        
        if(i==1)
            ylabel('Nuc Loc')
        end
  
        %concatenate all cells to calculate global max and min values (for
        %heatmap colorbar
        allCellsHeatMap = [allCellsHeatMap; varSync];
        indexes(idxAll: idxAll+size(varSync,1)-1) = i; %group vector
        idxAll= idxAll+size(varSync,1); 
         
    end
    
    minVal = min(min(allCellsHeatMap)); 
    maxVal = max(max(allCellsHeatMap)); 
    %Plot Heatmaps
    figure(1)
    for i=1:length(strains)
        subplot(nRows,length(strains),i+length(strains))
        heatmap(allCellsHeatMap(indexes==i,:),[],[],[],'colormap','parula','MinColorValue',minVal,'MaxColorValue',maxVal); 
        if(i==1)
            ylabel('Cell number')
        end
            xlabel('time')
        %number of cells in the title
        title(strcat('n=',num2str( sum(indexes==i))))
    end
  
    set(gcf,'color','w');
    set(gcf, 'PaperPositionMode', 'auto');
    suplabel(TFname,'t')

end

