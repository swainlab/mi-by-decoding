function varNorm = normMeanOffset(var,switchPoint,seriesLength,type)
%apply normalization function x_0(t)=( x(t)-mu )/sigma
%where x_0(t) is the normalised time series for a given cell
%and mu, sigma are mean and standard deviation for the
%population. Any other normalisation can be applied.

varmean = mean(mean(var(:,switchPoint-seriesLength-1:switchPoint-2)));
varsd = mean(std(var(:,switchPoint-seriesLength-1:switchPoint-2)));

varNorm = (var - varmean)/varsd;