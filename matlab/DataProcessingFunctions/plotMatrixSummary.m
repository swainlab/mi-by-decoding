function [ handles ] = plotMatrixSummary( errorMatrix ,xlab,ylab)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
colors = distinguishable_colors(size(errorMatrix,1)); 
hold on 
 boundedline(1:size(errorMatrix,2),mean(errorMatrix),std(errorMatrix),'black','alpha');
for i = 1:size(errorMatrix,1)
     handles(i) = plot(errorMatrix(i,:),'color',colors(i,:),'LineWidth',1)   ;
end

handles(i+1) = plot(mean(errorMatrix),'k','LineWidth',3);
ylabel(ylab)
xlabel(xlab); 
xlim([1,size(errorMatrix,2)])

end

