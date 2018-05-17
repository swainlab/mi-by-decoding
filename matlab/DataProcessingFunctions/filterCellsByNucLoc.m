function [nucLocFil,sumFil] = filterCellsByNucLoc(nucLoc,sdN,filterType,propTime,switchPoint)
%FILTERCELLSBYNUCLOC Filter out cells that display low signal
%nucLoc: trajectories of nuclear localisation
%sdN: standard deviation threshold for filering low-signal cells.
%filterType: if more than one filtering criteria is applied
%propTime: time window over which the filter will be calculated
%switchPoint: time point of the transition, used to consider the pre-shock
%
%Author: Alejandro Granados
%Copyright Gasper Tkacik and Peter Swain 2018

if(nargin<5)
    propTime = 0.8;
end

if(nargin<4)
    filterType = 1; %Default filter type, add more filters if necessary
end


%How many time points to consider when the filtering parameteres are
%calculated. Responses for transcription factors considered here last for
%20-50minutes (8-20timepoints)
respDuration = 10; %timepoints 2.5minutes

seriesLength  = 20; %total length of trajectories.
type = 2; %Type of normalisation that will be applied to the data after filtering

%Interpolate nans (if any) in the data.
nucLoc = inpaint_nans(nucLoc);
%apply normalization function x_0(t)=( x(t)-mu )/sigma
%where x_0(t) is the normalised time series for a given cell
%and mu, sigma are mean and standard deviation for the
%population. Any other normalisation can be applied.
nucLoc = normMeanOffset(nucLoc,switchPoint,seriesLength,type);

%calcualte mean and standard deviation in the pre-shock
%environment.
meanBefore = mean(nucLoc(:,switchPoint-respDuration));
stdBefore = std(nucLoc(:,switchPoint-respDuration));

lengthX = 2;  %how long to explore the time series for filtering decision (x =times the length used for MI calculation)

%Different filtering criteria can be applied depending on the
%data.
if(filterType==1)
    % filter based on low nuc Loc (e.g.,for sfp1)
    
    threshold = meanBefore - sdN*stdBefore;
    %remove cells for which the nucLoc is below the nStd during
    %the pre-shock segment of the time series.
    %only cells for which the nucLoc is below nStd for at least
    %propTime time points will be removed. Strong indicationd of
    %dead cells or false positive in segmentation.
    
    fil =  sum(nucLoc(:,switchPoint-seriesLength*lengthX:switchPoint)'<threshold)<(seriesLength*lengthX+1)*propTime;
    
    %Any other filtering/curation that suits the data set should be
    %applied here.
end


nucLocFil = nucLoc(fil,:);
nucLocFil = normMeanOffset(nucLocFil,switchPoint,seriesLength,type);

%How many cells were accepted by the filter.
sumFil=num2str(sum(fil));
end
