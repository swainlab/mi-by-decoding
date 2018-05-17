function  [infoDataTime, infoMatrix,errorMatrix]  = calculateInfoList(infdata,maxTP,bts)
%CALCULATEINFOLIST Plots the MI as a function of response duration
%
%Calculation of mutual information between the time-series stored in infdata
%and a two-state environmental transition (e.g., rich media & stress)
%Prints the MI as a function of the duration of the response and the
%confusion matrix (see SI). The confusion matrix is only shown for the
%timepoint for which the MI is maximum. MI and confusion matrix are
%averages across bootstrap replicates.
%
%infdata: structure with all experiments
%nBins: maximum number of timepoints to consider (from 1 to nBins)
%bts: number of bootstrap pseudoreplicates
%
%Authors: Alejandro Granados
%Copyright Gasper Tkacik and Peter Swain 2018

TFlist = fieldnames(infdata);

errorMatrix = zeros(length(TFlist),maxTP);
infoMatrix =  zeros(length(TFlist),maxTP);

mincolor = 0; maxcolor = 0.5;

%MI information usually converges for 5 principal components
%although more can be considered for complex data sets.
defaultNPC = 5;
s = 1; %only one repeat per tf in the list


for tf=1:length(TFlist)  %for each experiments
    for tp =1:maxTP %for increasing durations of the responses
        
        %just take as many PC as n of time points and 5
        %otherwise
        if(tp<defaultNPC), nPC = tp; else nPC = defaultNPC; end
        
        initial =  infdata.(TFlist{tf}).rich(:,1:tp);
        transient = infdata.(TFlist{tf}).stress(:,1:tp);
        
        %Save MI and confusion matrix for all (bts)
        %bootstrap replicates across all principal
        %components considered (nPC)
        %dim(M_aux)=[n_states,n_states,bts,nPC]
        
        %MAIN CALL TO METHOD
        [Ipca_aux,M_aux]=MIdecoding({initial,transient},bts,nPC);
        
        dat = Ipca_aux;
        %save all data for all time points.
        infoDataTime.(TFlist{tf}).(strcat('n',num2str(s))).timeData(tp).dat = dat;
        
        infoDataTime.(TFlist{tf}).M_time.timeData(tp).allM= M_aux; %
    end
end


for tf =1:length(TFlist)
    strain = TFlist{tf};
    infoTP = []; pcMax = []; infStd = []; decodingError = []; decodingErrorSTD=[];
    for tp =1:maxTP
        dat =infoDataTime.(strain).n1.timeData(tp).dat;
        [infoTP(tp),b] = max(nanmean(dat));
        pcMax(tp) = b; % PC with max info for this timepoint
        infStd(tp) = nanstd(dat(:,b));
        %Number of bootstraps
        nBoot = size(infoDataTime.(strain).M_time.timeData(tp).allM(:,:,:,b),3);
        
        %Calculate de decoding error
        e_ =zeros(1,nBoot);
        for nb = 1:nBoot
            %select the NPC for which information is max and then
            %calculate the mean error across bootstraps.
            thisM = infoDataTime.(strain).M_time.timeData(tp).allM(:,:,nb,b);
            e_(nb) = sum(sum((1-eye(2,2)).*thisM));
        end
        
        decodingError(tp)= mean(e_);  %average decoding error
        decodingErrorSTD(tp) = std(e_);
    end
    subplot(3,length(TFlist),tf);
    
    
    
    errorbar(infoTP,infStd); ylim([0,1]); title(strain); xlabel('time point');
    if(tf==1), ylabel('info (bits'); end
    
    [~,b] = max(infoTP);
    M = infoDataTime.(TFlist{tf}).M_time.timeData(b).allM;
    %Calculate the mean confusion matrix (across bootstraps) for the PC = pcMax(b)
    M_ = mean(M(:,:,:,pcMax(b)),3);
    
    subplot(3,length(TFlist),tf+length(TFlist));
    mlabel = arrayfun(@(x){sprintf('%0.3f',x)}, M_);   levelTag ={'u0','u1'};
    heatmap(M_, levelTag, levelTag, mlabel,'TickAngle', 45,'Colormap',@cool,'MinColorValue', mincolor, 'MaxColorValue', maxcolor);
    
    
    subplot(3,length(TFlist),tf+length(TFlist)*2)
    
    decodingError(1) = 0.5; decodingErrorSTD(1) = 0;
    
    errorbar(decodingError,decodingErrorSTD)
    ylim([0,1])
    %plot conf matrix con time point with max info
    %export data:
    infoMatrix(tf,:) = infoTP ;
    errorMatrix(tf,:) = decodingError;
    
end %end for TFs
end

