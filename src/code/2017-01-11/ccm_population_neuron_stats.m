function Data = ccm_population_neuron_stats(subject,projectRoot,projectDate,opt)

%
% function Data = ccm_population_neuron(subjectID, sessionID, plotFlag,
% unitArray)
%
% Population neuron collection for choice countermanding task. Only plots the
% sdfs. To see rasters, use ccm_single_neuron_rasters, which displays all
% conditions in a given epoch
%
% input:
%   subjectID: e.g. 'Broca', 'Xena', 'pm', etc
%   sessionID: e.g. 'bp111n01', 'Allsaccade'
%
%   opt: A structure with various ways to select/organize data: If
%   ccm_population_neuron.m is called without input arguments, the default
%   opt structure is returned. opt has the following fields with
%   possible values (default listed first):
%
%    opt.dataType = 'neuron', 'lfp', 'erp';
%
%    opt.doStops        = true, false;
%    opt.filterData 	= false, true;
%    opt.stopHz         = 50, <any number, above which signal is filtered;
%    opt.unitArray      = {'spikeUnit17a'},'each', units want to analyze
%
%
% Returns Data structure with fields:
%
%   Data.signalStrength(x).(condition).ssd(x).(epoch name)
%
%   condition can be:  goTarg, goDist, stopTarg, stopDist, stopStop
%   ssd(x):  only applies for stop trials, else the field is absent
%   epoch name: fixOn, targOn, checkerOn, etc.





%%
if nargin < 4
    opt.categoryName   	= 'presacc';
    opt.epochArray      = {'targOn','checkerOn','stopSignalOn','responseOnset'};
    opt.doStops         = true;
    opt.excludeSessions = true;
    opt.dataType        = 'neuron';
    opt.multiUnit    	= false;
    opt.normalize    	= false;
    opt.ms2Std          = 75;
    opt.ssrtUse         = 'intWeightPerSession';
    opt.saccadeBaseRatio         = [];
    if nargin == 0, Data = opt; return, end
end

if ~opt.doStops
    opt.epochArray      = {'targOn','checkerOn','responseOnset'};
end

if opt.multiUnit
    addMulti = '_multiUnit';
else
    addMulti = [];
end
if opt.normalize
    addNorm = '_normalized';
else
    addNorm = [];
end


% ____________________ CONSTANTS AND VARIABLES ____________________
% stopHz          = opt.stopHz;
% Define axes limits

goOutcomeArray      = {'goTarg'}; %{'goTarg', 'goDist'};
stopOutcomeArray    = {'stopTarg'}; %{'stopTarg', 'stopDist'};
conditionArray       = {'hardIn', 'hardOut', 'easyIn', 'easyOut'};
conditionArrayInd   = [1 2 3 4];

% ________________________________________________
% CHOICE AND COHERENCE (DDM-LIKE) CONSTNANTS AND PREPARE DATA

% These will determine the trial-to-trial epochs used for analyses:
epochOffset = 120;  % When to begin spike rate analysis after stimulus (checkerboard) onset
preSaccadeBuffer = 10; % When to cut off spike rate analysis before saccade onset
minEpochDuration = 10; % Only include trials for which the determined epoch is this long

epochRangeChecker = 1 : 1000;  % Make this big enough to catch really late cancel times.

% Initialize choice and coherence dependence to null assumption (i.e. they
% aren't)
choiceDependent     = false;
coherenceDependent  = false;
coherenceDependentRank  = false;
ddmLike             = false;
ddmLikeRank             = false;

% some constants
alphaChoice     = .05;   % alpha criteria for choice dependence
alphaCoherence  = .05;   % alpha criteria for coherence dependence





% ____________________    LOAD DATA    ____________________
dataPath = fullfile(projectRoot,'data',projectDate,subject);
load(fullfile(dataPath, ['ccm_',opt.categoryName,'_neurons',addMulti])) % Load neurons table for specific category population
if opt.excludeSessions
    sessionRemove = ccm_exclude_sessions(subject);
    neurons = neurons(~ismember(neurons.sessionID, sessionRemove),:);
end
if ~isempty(opt.saccadeBaseRatio)
    load(fullfile(dataPath, ['ccm_neuronTypes', addMulti]));
    includeInd = neuronTypes.saccadeBaseRatio >= opt.saccadeBaseRatio;
    keepUnit = neuronTypes(includeInd, 1:4);
    neurons = intersect(keepUnit, neurons);
end


Data = load(fullfile(dataPath, ['ccm_all_neuron_population',addMulti,addNorm])); % Load Data struct for whole population

% Index the relelvant units for the Go Data
goTargInd = ismember(cell2table(Data.unitGo,'VariableNames',{'sessionID' 'unit'}), neurons(:,1:2));



% load the population of cancel time anlysis
fileName = fullfile(dataPath, 'go_vs_canceled', opt.ssrtUse, ['ccm_canceled_vs_go_neuronTypes', addMulti]);
load(fileName)









% Get the presaccade : fixation baseline slope to establish criteria for
% inclusion
baselineEpoch = -99:0;
baselineRate = nanmean(nanmean(Data.targOn.easyIn.goTarg.sdf(goTargInd, Data.targOn.alignGo + baselineEpoch), 1));
saccadeEpoch = -50:0;
saccadeRate = nanmean(Data.responseOnset.easyIn.goTarg.sdf(goTargInd, Data.responseOnset.alignGo));

fprintf('%s\tSaccade / Baseline ratio: %.2f\n', opt.categoryName, saccadeRate / baselineRate);






%    CHOICE DEPENDENCE
% =================================================================


% Determine Choice Selection Time Method 1 -- via Hanes et al 1998 (p.822) differential sdf test
% ------------------------------
% Figure out when signals diverge w.r.t. choice direction. Use the
% easiest left and right coherence, so the epoch begins as early as
% possbile (to capture as much coherence divergence during
% coherence-dependence analysis
choiceSelectionTime = nan; % Initialize to NaN

goTargEasyInFn = nanmean(Data.checkerOn.easyIn.goTarg.sdf(goTargInd,:), 1);
goTargEasyOutFn = nanmean(Data.checkerOn.easyOut.goTarg.sdf(goTargInd,:), 1);

inOutDiffFn = goTargEasyInFn - goTargEasyOutFn;

% Get the standard deviation of the differential sdf during the
% 500 ms before checkerboard onset
stdDiff = std(inOutDiffFn(1 : Data.checkerOn.alignGo));

% GoIn-GoOut SDF from checker Onset to end of checker epoch.
sdfDiff = abs(inOutDiffFn(Data.checkerOn.alignGo : end))';


% are there times at which the difference between sdfs is
% greater than 2 standard deviations of the difference 500
% ms before checkerboard onset? Check from checkerboard onset
% to end of the checkerboard epoch.
% So std2Ind starts at checkerboard onset.
std2Ind = sdfDiff > 2*stdDiff;

% Look for a sequence of opt.ms2Std ms for which the go sdf is 2
% std greater than the stop sdf. Determein whether there was a time
% after the checkerboard onset that the differential
% sdf was > 2*Std for at least opt.ms2Std ms.
riseAbove2Std = find([0; diff(std2Ind)] == 1);
sinkBelow2Std = find([0; diff(std2Ind)] == -1);
if ~isempty(riseAbove2Std)
    % Get rid of occasions for which the signals differ
    % going into the epoch (and therefore they will
    % cease to differ before they begin again to
    % differ)
    if ~isempty(sinkBelow2Std)
        sinkBelow2Std(sinkBelow2Std < riseAbove2Std(1)) = [];
    end
    
    % If there's one more riseAbove2Std than sinkBelow2Std, the last riseAbove2Std
    % will last until the end of the sdf: Add to
    % singkBelowStd the end of the epoch
    if length(riseAbove2Std) > length(sinkBelow2Std)
        sinkBelow2Std = [sinkBelow2Std; length(goTargEasyInFn)-Data.checkerOn.alignGo];
    end
    
    % Now riseAbove2Std length should be equal. See if
    % any of the riseAbove2Std streaks go longer than
    % 50ms
    ind = find(sinkBelow2Std - riseAbove2Std >= opt.ms2Std, 1);
    if ~isempty(ind)
        choiceSelectionTime = riseAbove2Std(ind);
    end
end




% ________________________________________________________________
% DETERMINE THE ONSET OF THE CHOICE DEPENDENCE Method 2 a la Ding & Gold
tChoice = nan;  % Initialize to nan; in the case of eeg signals there may not be a tChoice (noise in signal)
nanThreshold = .3; % Establish a threshold of non-signal fraction over which we cease looking for tChoice
slideWindowWidth = 50; % ms, Ding and Gold used a 100ms sliding window
slideWindowStep = 2; % ms, D&G used 10 ms steps


choiceDependenceFound = false;
iStepInd = 0;
while ~choiceDependenceFound
    
    
    iEpochBegin = (Data.checkerOn.alignGo + (iStepInd * slideWindowStep) + epochOffset) * ones(sum(goTargInd), 1);
    iEpochEnd = iEpochBegin + slideWindowWidth;
    
    
    if iEpochEnd(1) > length(goTargEasyInFn)
        break
    end
    % Choice dependence
    inMetric = mean(Data.checkerOn.easyIn.goTarg.sdf(goTargInd,iEpochBegin:iEpochEnd), 2);
    outMetric = mean(Data.checkerOn.easyOut.goTarg.sdf(goTargInd,iEpochBegin:iEpochEnd), 2);
    
    [p, h, stats] = ranksum(inMetric, outMetric);
    
    if p < alphaChoice
        choiceDependenceFound = true;
        tChoice = (iStepInd * slideWindowStep) + epochOffset;  % currently use choiceSelectionTime instead of tChoice
    end
    iStepInd = iStepInd + 1;
end

choiceSelectionTime = tChoice;




% Define the beginning and end of the epoch to analyze for choice and coherence
% dependence
% Initialize end of epoch as RTs from easiect choice
% conditions
epochEnd = ceil(Data.checkerOn.easyIn.goTarg.rt(goTargInd));
% Use Choice selection time if possible, to define
% beginning of epoch. Otherwise, use the RTs
if ~isnan(choiceSelectionTime)
    epochBegin = choiceSelectionTime * ones(sum(goTargInd), 1);
    %     epochBegin = 175 * ones(sum(goTargInd), 1);
else
    epochBegin = ceil(.5 * epochEnd);
end
% epochEnd = 100 + epochBegin;


% Adjust for the alignment index
epochEnd = epochEnd + Data.checkerOn.alignGo;
epochBegin = epochBegin + Data.checkerOn.alignGo;
epochDuration = epochEnd - epochBegin;

% If there are units with negative epochs because of the choice selection
% time w.r.t RT, adjust them
negativeEpochUnit = epochEnd < epochBegin + minEpochDuration;
epochBegin(negativeEpochUnit) = epochEnd(negativeEpochUnit) - minEpochDuration;


% sums of the average rasters per unit:
inMetric = sum(Data.checkerOn.easyIn.goTarg.raster(goTargInd,epochBegin:epochEnd), 2);
outMetric = sum(Data.checkerOn.easyOut.goTarg.raster(goTargInd,epochBegin:epochEnd), 2);
% inMetric = mean(Data.checkerOn.easyIn.goTarg.sdf(goTargInd,epochBegin:epochEnd), 2);
% outMetric = mean(Data.checkerOn.easyOut.goTarg.sdf(goTargInd,epochBegin:epochEnd), 2);



% Choice dependence

[p, h, stats]   = signrank(inMetric, outMetric);

if p < alphaChoice
    choiceDependent = true;
end






%    COHERENCE DEPENDENCE
% =================================================================

% For IN trials
% inSpikeRate     = [mean(Data.checkerOn.hardIn.goTarg.sdf(goTargInd,epochBegin:epochEnd), 2); mean(Data.checkerOn.easyIn.goTarg.sdf(goTargInd,epochBegin:epochEnd), 2)];
inSpikeRate     = [sum(Data.checkerOn.hardIn.goTarg.raster(goTargInd,epochBegin:epochEnd), 2); sum(Data.checkerOn.easyIn.goTarg.raster(goTargInd,epochBegin:epochEnd), 2)];
inCond          = [ones(sum(goTargInd), 1); 2 * ones(sum(goTargInd), 1)];
% For OUT trials
% outSpikeRate     = [mean(Data.checkerOn.hardOut.goTarg.sdf(goTargInd,epochBegin:epochEnd), 2); mean(Data.checkerOn.easyOut.goTarg.sdf(goTargInd,epochBegin:epochEnd), 2)];
outSpikeRate     = [sum(Data.checkerOn.hardOut.goTarg.raster(goTargInd,epochBegin:epochEnd), 2); sum(Data.checkerOn.easyOut.goTarg.raster(goTargInd,epochBegin:epochEnd), 2)];
outCond          = [ones(sum(goTargInd), 1); 2 * ones(sum(goTargInd), 1)];



% Regress spikeRate vs signalStrength into RF
[coeffIn, sIn]          = polyfit(inCond, inSpikeRate, 1);
[yPredIn, deltaIn]  = polyval(coeffIn, inCond, sIn);
statsIn             = regstats(inCond, inSpikeRate);
rIn                 = corr(inCond, inSpikeRate);

slopeIn     = coeffIn(1);
signSlopeIn = sign(slopeIn);
fTestIn     = statsIn.fstat.f;
pValIn      = statsIn.fstat.pval;

% Regress spikeRate vs signalStrength out of RF
[coeffOut, sOut]        = polyfit(outCond, outSpikeRate, 1);
[yPredOut, deltaOut] = polyval(coeffOut, outCond, sOut);
statsOut            = regstats(outCond, outSpikeRate);
rOut                = corr(outCond, outSpikeRate);

slopeOut    = coeffOut(1);
signSlopeOut = sign(slopeOut);
fTestOut    = statsOut.fstat.f;
pValOut     = statsOut.fstat.pval;




% Decision tree to determine whether the neuron/signal was "coherence dependent"
if pValIn < alphaCoherence && pValOut > alphaCoherence
    coherenceDependent = true;
elseif pValIn > alphaCoherence && pValOut < alphaCoherence
    coherenceDependent = true;
elseif pValIn < alphaCoherence && pValOut < alphaCoherence
    % slopeOut must have opposite sign than slopeIn
    if signSlopeIn ~= signSlopeOut
        coherenceDependent = true;
    end
end

% if choiceDependent && coherenceDependent
%     ddmLike = true;
% end






% Alternate Coherence Dependence test- Rank sum, like choice dependence
% In RF Easy vs Hard
% [pIn, h, stats]   = ranksum(mean(Data.checkerOn.easyIn.goTarg.sdf(goTargInd,epochBegin:epochEnd), 2), mean(Data.checkerOn.hardIn.goTarg.sdf(goTargInd,epochBegin:epochEnd), 2));
% [pIn, h, stats]   = ranksum(sum(Data.checkerOn.easyIn.goTarg.raster(goTargInd,epochBegin:epochEnd), 2), sum(Data.checkerOn.hardIn.goTarg.raster(goTargInd,epochBegin:epochEnd), 2));
% [h,p,ci,stats] = ttest(sum(Data.checkerOn.easyIn.goTarg.raster(goTargInd,epochBegin:epochEnd), 2), sum(Data.checkerOn.hardIn.goTarg.raster(goTargInd,epochBegin:epochEnd), 2))% Out of RF Easy vs Hard
[pIn,h,stats]   = signrank(sum(Data.checkerOn.easyIn.goTarg.raster(goTargInd,epochBegin:epochEnd), 2), sum(Data.checkerOn.hardIn.goTarg.raster(goTargInd,epochBegin:epochEnd), 2));
% [pOut, h, stats]   = ranksum(mean(Data.checkerOn.easyOut.goTarg.sdf(goTargInd,epochBegin:epochEnd), 2), mean(Data.checkerOn.hardOut.goTarg.sdf(goTargInd,epochBegin:epochEnd), 2));
% [pOut, h, stats]   = ranksum(sum(Data.checkerOn.easyOut.goTarg.raster(goTargInd,epochBegin:epochEnd), 2), sum(Data.checkerOn.hardOut.goTarg.raster(goTargInd,epochBegin:epochEnd), 2));
[pOut, h, stats]   = signrank(sum(Data.checkerOn.easyOut.goTarg.raster(goTargInd,epochBegin:epochEnd), 2), sum(Data.checkerOn.hardOut.goTarg.raster(goTargInd,epochBegin:epochEnd), 2));

% Decision tree to determine whether the neuron/signal was "coherence dependent"
if pIn <= alphaCoherence && pOut > alphaCoherence
    coherenceDependentRank = true;
elseif pIn > alphaCoherence && pOut <= alphaCoherence
    coherenceDependentRank = true;
elseif pIn <= alphaCoherence && pOut <= alphaCoherence
    % slopeOut must have opposite sign than slopeIn
    if signSlopeIn ~= signSlopeOut
        coherenceDependentRank = true;
    end
end

if choiceDependent && coherenceDependentRank
    ddmLike = true;
end












%  STOPPING/CANCELING POPULATION STATS


% Index the relelvant units for the StopStop Data
unitIndEasy = ismember(cell2table(Data.checkerOn.easyIn.stopStop.unitStop,'VariableNames',{'sessionID' 'unit'}), neurons(:,1:2));
unitIndHard = ismember(cell2table(Data.checkerOn.hardIn.stopStop.unitStop,'VariableNames',{'sessionID' 'unit'}), neurons(:,1:2));
% Concatenate all data from the same units
cancelEasy = cell2table(Data.checkerOn.easyIn.stopStop.unitStop(unitIndEasy,:),'VariableNames',{'sessionID' 'unit'});
cancelEasy.ssd = Data.checkerOn.easyIn.stopStop.ssd(unitIndEasy);
cancelHard = cell2table(Data.checkerOn.hardIn.stopStop.unitStop(unitIndHard,:),'VariableNames',{'sessionID' 'unit'});
cancelHard.ssd = Data.checkerOn.hardIn.stopStop.ssd(unitIndHard);


% Build a table of ssrt and cancel time
cancelEasy.ssrt = nan(length(cancelEasy.ssd), 1);
cancelEasy.cancelTime = nan(length(cancelEasy.ssd), 1);
cancelHard.ssrt = nan(length(cancelHard.ssd), 1);
cancelHard.cancelTime = nan(length(cancelHard.ssd), 1);

for i = 1 : length(cancelEasy.ssd)
    iUnitInd = ismember(cancelTypes.sessionID, cancelEasy.sessionID(i)) & ismember(cancelTypes.unit, cancelEasy.unit(i));
    iSsdInd = cancelTypes.stopStopSsd{iUnitInd} == cancelEasy.ssd(i);
    iCondInd = cancelTypes.stopStopCond{iUnitInd} == 1;
    cancelEasy.ssrt(i) = cancelTypes.stopStopSsrt{iUnitInd}(iSsdInd & iCondInd);
    cancelEasy.cancelTime(i) = cancelTypes.cancelTime2Std{iUnitInd}(iSsdInd & iCondInd) - cancelEasy.ssrt(i) - cancelEasy.ssd(i);
end
for i = 1 : length(cancelHard.ssd)
    iUnitInd = ismember(cancelTypes.sessionID, cancelHard.sessionID(i)) & ismember(cancelTypes.unit, cancelHard.unit(i));
    iSsdInd = cancelTypes.stopStopSsd{iUnitInd} == cancelHard.ssd(i);
    iCondInd = cancelTypes.stopStopCond{iUnitInd} == 2;
    cancelHard.ssrt(i) = cancelTypes.stopStopSsrt{iUnitInd}(iSsdInd & iCondInd);
    cancelHard.cancelTime(i) = cancelTypes.cancelTime2Std{iUnitInd}(iSsdInd & iCondInd) - cancelHard.ssrt(i) - cancelHard.ssd(i);
end















return


