function [MutInf,confM,Nerrs] = MIdecoding(tcdata,varargin)
%MIDECODING Estimate mutual information between states and time series
%   MutInf = MIDECODING(tcdata) Calculates a confusion matrix based on the 
%   performance of a classifier function and uses this to estimate a lower
%   bound on the mutual information (MI) between states and time series.
%   - tcdata: a struct or cell array with fields/elements containing time
%       series for each state (cell index is dim 1; time index is dim 2).
%       Note that all time courses must have equal numbers of time points, 
%       but can differ in the number of cells.
%   - MutInf: a matrix of MI estimates with bootstraps along dim 1 and an
%       increasing number of included components (or time points if 'PCA'
%       is set to 'none') along dim 2.
%
%   [MutInf,confM,Nerrs] = MIDECODING(tcdata) returns also:
%   - confM: confusion matrices for a classification of S states obtained 
%       for N bootstraps and up to M included components in an 
%       S x S x N x M array.
%   - Nerrs: total number of errors made by the classifier in an N x M
%       matrix.
%
%   MIDECODING(tcdata,N_bootstrap,N_components) allows for optional
%   specification of the following parameters:
%   - N_bootstrap (default 25): number of bootstrap iterations to run
%   - N_components (default 10): maximum number of PCA components to
%       consider (or the maximum number of time points if 'PCA' is set to
%       'none')
%   
%   MIDECODING(...,'classifier',classifier) 
%       Specify the classifier to use. Defaults to 'svmlin-onevsone', a
%       linear SVM classifier using the one-vs-one method for multi-class
%       classification, but can be changed to 'randomforest' for a 
%       random forest classifier (using Matlab's TreeBagger class), to 
%       'svmlin-complete' for complete multi-class linear SVM 
%       classification, or to 'svmlin-dendrogram' for a multi-class
%       dendrogram-based linear SVM. Note that if tcdata contains only two 
%       states, then the SVM methods all revert to the standard two-class 
%       SVM.
%
%   MIDECODING(...,'PCA',type) Specify the type of PCA to perform. 
%       'train' (the default) uses the training data to perform PCA, and 
%       'none' turns off the PCA.
%   
%   MIDECODING(...,'min_components',m) Specify the minimum number of PCA
%       components to include.
%
%   MIDECODING(...,'normalise',type) Specify the type of normalisation to
%       perform. 'none' (the default) performs no normalisation (except for
%       that automatically done by the classifier), 'raw' performs
%       normalisation based on the complete data set, and 'bootstrap'
%       performs normalisation on the reduced bootstrap data sets.
%
%   MIDECODING(...,'softnorm',false) If normalisation is applied, then use
%       the maximum and minimum to perform a 'hard' normalisation rather
%       than the standard deviation and mean of the default 'soft'
%       normalisation.
%
%   MIDECODING(...,'featnorm',false) If normalisation is applied, then
%       apply it uniformly to all time points (i.e., the features), rather 
%       than the default of normalising each time point independently.
%
%   MIDECODING(...,'cost',cost) Set the cost parameter used by the SVM
%       classifiers.
%   
%   If you publish results that make use this software or the mutual 
%   information by decoding algorithm, please cite:
%   Granados, A.A., Pietsch, J.M.J., Cepeda-Humerez, S.A., Farquhar, I.L.,
%   Tkacik, G., and Swain, P.S. (2018) Distributed and dynamic
%   intracellular organization of extracellular information.
%
%   Authors: Sarah Cepeda, Alejandro Granados and Julian Pietsch
%   Copyright Gasper Tkacik and Peter Swain 2018

if isstruct(tcdata) && isscalar(tcdata)
    tcdata = struct2cell(tcdata);
end

ip = inputParser;
ip.addRequired('tcdata',@(x) iscell(x) && isvector(x) && all(cellfun(@isnumeric,x)));
ip.addOptional('N_bootstrap',25,@(x) isnumeric(x) && isscalar(x));
ip.addOptional('N_components',10,@(x) isnumeric(x) && isscalar(x));
ip.addParameter('classifier','svmlin-onevsone',@(x) ischar(x));
ip.addParameter('min_components',1,@(x) isnumeric(x) && isscalar(x));
ip.addParameter('uniform',true,@(x) isscalar(x) && islogical(x));
ip.addParameter('normalise','none',@(x) ischar(x));
ip.addParameter('softnorm',true,@(x) isscalar(x) && islogical(x));
ip.addParameter('featurenorm',true,@(x) isscalar(x) && islogical(x));
ip.addParameter('PCA','train',@(x) ischar(x) && ismember(x,{'train','none'}));
ip.addParameter('cost',1,@(x) isscalar(x) && isnumeric(x));
ip.parse(tcdata,varargin{:});

N_bootstrap = ip.Results.N_bootstrap;
N_components = ip.Results.N_components;
min_components = ip.Results.min_components;
classifier = ip.Results.classifier;
uniform = ip.Results.uniform;
normalise = ip.Results.normalise;
softnorm = ip.Results.softnorm;
featurenorm = ip.Results.featurenorm;
PCAtype = ip.Results.PCA;
cost = ip.Results.cost;

% Raise warning if NaNs encountered in data:
if any(cellfun(@(x) any(isnan(x(:))),tcdata))
    warning('NaNs found in data for MI calculation. Unexpected results may arise.');
end

% Raise error if the data arrays have conflicting numbers of time points:
N_timepoints = cellfun(@(x) size(x,2),tcdata);
if all(N_timepoints==N_timepoints(1))
    N_timepoints = N_timepoints(1);
else
    error('All classes must have the same number of time points');
end

% Always limit the number of components to the number of time points
N_components = min(N_timepoints,N_components);

N_classes = length(tcdata);
if N_classes==2 && ~strcmp(classifier,'randomforest')
    classifier = 'svmlin-2state';
end

% Use only 90% of the available cells for each bootstrap
N_cells = cellfun(@(c) floor(0.9*size(c,1)),tcdata);
if uniform, N_cells = repmat(min(N_cells),size(N_cells)); end

% Use only a quarter of the bootstrap sample for testing
N_test = round(1/4*N_cells);
N_train = N_cells - N_test;

if strcmp(normalise,'raw')
    tcdata = calcnorm(tcdata,softnorm,featurenorm);
end

confM = NaN(N_classes,N_classes,N_bootstrap,N_components); % Confusion matrix
MutInf = NaN(N_bootstrap,N_components); % Information
Nerrs = NaN(N_bootstrap,N_components); % Number of errors in each trial

%for i = 1:N_bootstrap
parfor i = 1:N_bootstrap
    
    if strcmp(classifier,'svmlin-dendrogram')
        % Disable warning raised by pdist during calculation of the dendrogram:
        warning('off','stats:pdist:ConstantPoints');
        % Note that this warning occurs when only one vector component is used
        % to create the tree. It needs to be called within the parfor loop
        % since it must be set for each parallel process separately.
    end
    
    % Sample cells from each time course for this bootstrap
    data_i = cellfun(@(c,n) c(randperm(size(c,1),n),:),...
        tcdata,num2cell(N_cells),'Uniform',false);
    
    % Split data into training and test sets
    train_data_i = cellfun(@(d,n) d(1:n,:),...
        data_i,num2cell(N_train),'Uniform',false);
    test_data_i = cellfun(@(d,nte,ntr) d(ntr+1:ntr+nte,:),...
        data_i,num2cell(N_test),num2cell(N_train),'Uniform',false);
    
    pcaTF = [];
    if ~strcmp(PCAtype,'none')
        % Calculate principal components for this sample
        if strcmp(PCAtype,'train')
            pcaTF = pca(vertcat(train_data_i{:}));
        else
            error('Unrecognised PCA method encountered.');
        end

        % Project both test and training data onto principal components
        train_data_i = cellfun(@(c) c*pcaTF,train_data_i,'Uniform',false);
        test_data_i = cellfun(@(c) c*pcaTF,test_data_i,'Uniform',false);
    end
    
    if strcmp(normalise,'bootstrap')
        % Determine normalisation on the training data, and scale test_data
        % in the same way:
        [train_data_i,test_data_i] = calcnorm(train_data_i,...
            softnorm,featurenorm,test_data_i);
    end
    
    % Pre-initialise M just in case it gets used before it is set by the
    % try-catch loop (also avoids warning about temporary variables in
    % parallel for loops):
    Mij = NaN(N_classes,N_classes);
    
    % Loop over PCA components
    for j=min_components:N_components
        if strcmp(classifier,'svmlin-dendrogram')
            % NB: for complex hierarchies the dendrogram classifier can
            % fail sometimes, so catch these errors and leave as NaNs
            try
                [Mij,Nerrs(i,j)] = MultiSVMconfM(train_data_i,test_data_i,j,classifier,cost);
            catch
                continue
            end
        else
            [Mij,Nerrs(i,j)] = MultiSVMconfM(train_data_i,test_data_i,j,classifier,cost);
        end
        MutInf(i,j) = Info(Mij);
        confM(:,:,i,j) =  Mij;
    end
end

end

function [M,e] = MultiSVMconfM(train_data,test_data,k,method,cost)
N_classes = length(train_data);
if length(train_data)~=length(test_data)
    error('Training and test data sets must have equal numbers of classes');
end

% Generate a class label for each cell
classlabels = 0:N_classes-1;
N_train = cellfun(@(d) size(d,1),train_data);
N_test = cellfun(@(d) size(d,1),test_data);
train_label = arrayfun(@(l,n) l*ones(n,1),...
    classlabels(:),N_train(:),'Uniform',false);
test_label = arrayfun(@(l,n) l*ones(n,1),...
    classlabels(:),N_test(:),'Uniform',false);

if ~strcmp(method,'randomforest')
    % Subset the data to the specified number of components:
    train_data = cellfun(@(d) d(:,1:k),train_data,'Uniform',false);
    test_data = cellfun(@(d) d(:,1:k),test_data,'Uniform',false);
end

% Concatenate the test data into a single matrix
test_mat = vertcat(test_data{:});
test_label = vertcat(test_label{:});

switch method
    case 'svmlin-2state'
        MaxIter = 5000000;
        svm_obj = fitcsvm(vertcat(train_data{:}),vertcat(train_label{:}),...
            'IterationLimit',MaxIter,'KernelFunction','linear',...
            'Cost',cost*(1-eye(2)),'Standardize',true);
        class_test = svm_obj.predict(test_mat);
    case 'randomforest'
        bagger = TreeBagger(k*10,vertcat(train_data{:}),vertcat(train_label{:}));
        class_test = cellfun(@str2double,bagger.predict(test_mat));
    case 'svmlin-complete'
        ecoc = fitcecoc(vertcat(train_data{:}),vertcat(train_label{:}),...
            'Coding','binarycomplete','Cost',cost*(1-eye(N_classes)),...
            'Learners',templateSVM('Standardize',true));
        class_test = ecoc.predict(test_mat);
    case 'svmlin-onevsone'
        ecoc = fitcecoc(vertcat(train_data{:}),vertcat(train_label{:}),...
            'Coding','onevsone','Cost',cost*(1-eye(N_classes)),...
            'Learners',templateSVM('Standardize',true));
        class_test = ecoc.predict(test_mat);
    case 'svmlin-dendrogram'
        [svmstruct,level] = Train_DSVM(train_data',train_label','linear');
        if length(unique(level))~=length(level)
            error('MutualInformationCode:MultiSVM_PCA','Ambiguous levels...');
        end
        [class_test] = Classify_DSVM(test_mat,classlabels,svmstruct,level);
    otherwise
        error('Unknown classification method');
end

% Total number of errors in test
e = sum(class_test(:)~=test_label(:));

% Confusion matrix
M = zeros(N_classes);
for rl = 1:N_classes % loop over real class
    for est = 1:N_classes % loop over estimated class
        M(rl,est) = sum(classlabels(est) == ...
            class_test(test_label==classlabels(rl)));
    end
end
M = M/length(class_test);
end

function [normdata,normtestdata] = calcnorm(data,soft,perfeature,testdata)
if nargin<4
    testdata = {};
end
normtestdata = {};

% Normalisation stats must be calculated on the combined data
alldata = vertcat(data{:});
if soft
    % Normalise using mean as offset and 2x std. dev. to scale
    if perfeature
        % Normalise each column independently
        offset = mean(alldata,1);
        scale = 2*std(alldata,0,1);
    else
        % Normalise all columns by the same value
        offset = repmat(mean(alldata(:)),1,size(alldata,2));
        scale = repmat(std(alldata(:)),1,size(alldata,2));
    end
else
    % Normalise to [-1,1] using min and max
    if perfeature
        % Normalise each column independently
        minval = min(alldata,[],1);
        maxval = max(alldata,[],1);
    else
        % Normalise all columns by the same value
        minval = repmat(min(alldata(:)),1,size(alldata,2));
        maxval = repmat(max(alldata(:)),1,size(alldata,2));
    end
    offset = (maxval+minval)/2;
    scale = (maxval-minval)/2;
end
% If there are any columns with zero variance, set the scale to 1 to avoid
% division by zero:
scale(scale==0) = 1;
% Replicate data for all rows
normdata = cellfun(@(x) ...
    (x-offset(ones(size(x,1),1),:))./scale(ones(size(x,1),1),:),...
    data,'Uniform',false);
if ~isempty(testdata)
    normtestdata = cellfun(@(x) ...
        (x-offset(ones(size(x,1),1),:))./scale(ones(size(x,1),1),:),...
        testdata,'Uniform',false);
end
end