% Demonstration script for estimation of mutual information by decoding as 
% described in:
%   Granados, A.A., Pietsch, J.M.J., Cepeda-Humerez, S.A., Farquhar, I.L.,
%   Tkacik, G., and Swain, P.S. (2018) Distributed and dynamic
%   intracellular organization of extracellular information.
% In this script, mutual information between single-cell time series data 
% and environmental conditions is estimated.
%
% If you publish results that make use this software or the mutual 
% information by decoding algorithm, please cite the reference above.
% Experimental data from the paper can be downloaded from 
% http://dx.doi.org/10.7488/ds/2214
%
%   Authors: Alejandro Granados and Julian Pietsch
%   Copyright Gasper Tkacik and Peter Swain 2018

%% 2-state example of mutual information calculation
% Trajectories represent the nuclear localisation of Sfp1 in response to a
% transition from environment A (high glucose) to enviornment B (low
% glucose). The distribution of enviromental states is discrete and, for 
% the purpose of this example, assumed to be uniform (see main text and SI). 

% Load experimental data from figure 1. This step requires the JSONlab
% toolbox to be installed.
file_string = fileread('fig1_sfp1_replicates.json');
expts = loadjson(regexprep(file_string,'null','NaN'));

% Align the data to a common time and partition into environment A and B
infdata = struct();
reps = fieldnames(expts);


%parameters for filtering unhealthy/dead/noisy cells
filterLow = 1; nStd = 0.5; propTime=0.8;  

% Data pre-processing for each experiment
% In this example we analyze 6 experimental replicates of Sfp1 during a
% transition from high (2%) to low (0.1%) glucose.
for r=1:length(reps)
    times = expts.(reps{r}).general.times;
    midtimes = median(times);
    var = expts.(reps{r}).GFP.nucLoc;
    varSync = NaN(size(var));
    
    % Use linear interpolation to align traces to a common time
    for c=1:size(var,1)
        mask = ~isnan(var(c,:));
        varSync(c,:) = interp1(times(c,mask),var(c,mask),midtimes,'linear');                
    end
    
    % Time of the transition into stress
    origin = expts.(reps{r}).general.origin;
    
    % In rich media, Sfp1 is mostly localised in the nucleus albeit showing 
    % also stochastic translocation events. The filter in this example is 
    % specific for Sfp1
    
    % Filter cells with extreme low values of Sfp1 nuclear localisation in 
    % rich media (indication of unhealthy/dead cells or errors in 
    % segmentation), save the number of cells that remain in the dataset in
    % nCells
    [varSync,nCells{r}] = filterCellsByNucLoc(varSync,nStd,filterLow,propTime,origin); 
    
    %Select the regions in the data set that will be classified
    %In this experiment "rich" is the portion of the time series before the
    %stress was applied (at time=origin). Stress is the portion of the time
    %series after the stress was applied. 20 time points are considered
    %here.
    infdata.(reps{r}).rich = varSync(:,origin-20+(0:19));
    infdata.(reps{r}).stress = varSync(:,origin+(0:19));
end

%%
%plot(nanmean(expts.rep1.cy5.imBackground))
figure(2); clf; hold on;
plot(nanmean(expts.rep1.GFP.nucLoc(:,expts.rep1.general.origin-20+(0:19))))
plot(nanmean(expts.rep1.GFP.nucLoc(:,expts.rep1.general.origin+(0:19))))


%% DATA inspection
%visualize a summary of the data after pre-processing. 
%number of single cells to add in the plot. Data normalisation allows for
%comparing experimental replicates
nSingles=20; 
visualizeExper( infdata, [-20,20],[-3,3],nSingles,'Sfp1',fieldnames(infdata)) ;

%% Calculation of the 2-way MI
%Calculate the 2-way mutual information (MI) between the trajectories and the discrete
%environment (2 possible states: [stress,rich]; max_MI = 1bit)
%MI is calculated as function of the response duration (by considering increasing length of the trajectories)

%bts: Number of bootstrap pseudoreplicates (MI usually converges for bts>30 bootstrap pseudoreplicates)
bts = 30; 
%seriesLength:Maximum length of the time series
seriesLength=20;                                        
%infoMatrix & errorMatrix are the MI and classification error for each experiment as a function of the response's
%duration
figure,
[infoDataTime,infoMatrix,errorMatrix] =calculateInfoList(infdata,seriesLength,bts); 

%optional 
%Calculate the MI not as a function of the duration but considering the
%whole time series.
%mutinf = structfun(@(x) MIdecoding({x.rich,x.stress}),infdata,'Uniform',false);


%% Plot MI and decoding error
%Fraction of correctly classified (Fraction correct) and missclassified cells (decoding error)
%The equivalent MI is plotted for each experiment 
%(max_MI=1 bit correspondst to error=0; min_MI=0 bit corresponds to error=0.5, see SI) 
%The relation between decoding error and MI is explained in detail in the Supplemetary information

%Mean and standard deviation across experimental replicates is shown as shaded area (compare with Fig1 of the manuscript)
legendN = strcat('n=',nCells); 
legendN{end+1} = 'Mean'; 
figure,
subplot(2,2,2)
plotMatrixSummary(infoMatrix,'Response duration','2-way MI') ; 
ylim([0,1])
%plot mean information with individual experiments: 
% 
subplot(2,2,3) %plot error as function of duration. 
handles =plotMatrixSummary(errorMatrix,'Response duration','Error probability');
legend(handles,legendN)
ylim([0,1])
% prediction probability
subplot(2,2,1) %plot error as function of duration. 
handles =plotMatrixSummary(1-errorMatrix,'Prediction probability','Fraction correct');
ylim([0,1])


