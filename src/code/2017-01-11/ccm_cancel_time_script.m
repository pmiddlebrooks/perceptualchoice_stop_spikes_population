%%

subject = 'joule';
% subject = 'joule';

sessionRemove = ccm_exclude_sessions(subject);

saccadeBaseRatio = [];
saccadeBaseRatio = 2;

category = 'presacc_cancel_meanDifference';
category = 'presacc_ddmRankMeanStim_cancel_meanDifference';
% category = 'presacc_cancel_trialByTrial';
% category = 'presacc_ddmRankMeanStim_cancel_trialByTrial';
% category = 'presacc_cancel_meanSdf';
% category = 'presacc_ddmRankMeanStim_cancel_meanSdf';

nTrialCriteria = 10;
nTrialCriteria = 20;
nTrialCriteria = 1020;




projectDate = '2017-01-11';
projectRoot = '~/perceptualchoice_stop_spikes_population';

addpath(genpath(fullfile(projectRoot,'src/code',projectDate)));
dataPath = fullfile(projectRoot,'data',projectDate,subject);


% load the population of cancel time anlysis
fileName = fullfile(dataPath, 'go_vs_canceled', ssrtUse, ['ccm_canceled_vs_go_neuronTypes', addMulti]);
load(fileName)


fileName = fullfile(dataPath, ['ccm_',category,'_neurons', addMulti]);
load(fileName);

neurons = neurons(~ismember(neurons.sessionID, sessionRemove),:);
if ~isempty(saccadeBaseRatio)
    load(fullfile(dataPath, ['ccm_neuronTypes', addMulti]));
    includeInd = neuronTypes.saccadeBaseRatio >= saccadeBaseRatio;
    keepUnit = neuronTypes(includeInd, 1:4);
    neurons = intersect(keepUnit, neurons);
end


cancelData = table();
for j = 1 : size(neurons, 1)
    
    % find the indices in cancelTypes that correspond to this unit
    iInd = strcmp(neurons.sessionID(j), cancelTypes.sessionID) & strcmp(neurons.unit(j), cancelTypes.unit);
    find(iInd)
    if nTrialCriteria == 10
        cancelData = [cancelData; cancelTypes(iInd,:)];
    elseif nTrialCriteria == 20
        % Only use units that have conditions with at least 20 trials
        if sum(cancelTypes.nStopStop{iInd} >= 20)
            cancelData = [cancelData; cancelTypes(iInd,:)];
        end
    elseif nTrialCriteria == 1020
        if sum(cancelTypes.nStopStop{iInd} >= 20)
            j20Ind = cancelTypes.nStopStop{iInd} >= 20;
            for k = 4 : size(cancelTypes, 2)
                %             cancelTypes{iInd,k}
                % k
                %             cancelTypes{iInd,k}{:}(:) = cancelTypes{iInd,k}{:}(j20Ind);
                cancelTypes{iInd,k}{:}(~j20Ind) = [];
            end
        end
            cancelData = [cancelData; cancelTypes(iInd,:)];
    end
    
end

size(cancelData)
%%
% 20+ trials. Find the conditions with 10+ or 20+ trials
% validInd = cellfun(@(x) x >= nTrialCriteria, cancelData.nStopStop, 'uni', false);
% for k = 4 : size(cancelData, 2)
%     %     kCell = table2cell(cancelData(:,k));
%     cancelData(:,k) = cellfun(@(in1, in2) in1(in2), table2cell(cancelData(:,k)), validInd, 'uni', false);
% end

%% Old Cancel Time
cancelTime = cellfun(@(x,y,z) x - y - z, cancelData.cancelTime2Std, cancelData.stopStopSsd, cancelData.stopStopSsrt, 'uni', false);

% Grand average cancel time
% ssrtGrandAvg = mean(cell2mat(cancelData.stopStopSsrt));
% cancelTime = cellfun(@(x,y) x - y - ssrtGrandAvg, cancelData.cancelTime2Std, cancelData.stopStopSsd, 'uni', false);
% 
cancelData.cancelTime = cancelTime;
%% Keep only the conditions that cancel?
cancelLogical = cellfun(@(in1) in1 <= 20, cancelData.cancelTime, 'uni', false);
for j = 4 : size(cancelData, 2)
    % for k = 1 : size(cancelData, 1)
    cancelData(:,j) = cellfun(@(in1, in2) in1(in2), table2cell(cancelData(:,j)), cancelLogical, 'uni', false);
    % end
end
cancelTime = cancelData.cancelTime;

%% Trial by Trial sdf deflection Cancel Time


% % A neuron "cancels" if one of its condition's mean cancel time is within latestTime
% cancelTimeDistMean = cell(size(cancelData, 1), 1);
% for j = 1 : size(cancelData, 1)
%     cancelTimeDistMean{j} = cellfun(@nanmean, cancelData.cancelTimeDist{j});
% end
% cancelData.cancelTimeDistMean = cancelTimeDistMean;
%
% cancelTime = cell2mat(cancelData.cancelTimeDistMean);

%% Mean Sdf deflection cancel time
% cancelTime = cell2mat(cancelData.cancelTimeSdf);


%% Remove nans
% nStopStop = cell2mat(cancelData.nStopStop);
%
% nanCancel = isnan(cancelTime);
% nStopStop = nStopStop(~nanCancel);
% cancelTime = cancelTime(~nanCancel);

%% Weighted mean
weightCancelTime = cell2mat(cancelTime);
nStopStop = cell2mat(cancelData.nStopStop);
ridNan = isnan(weightCancelTime);
weightCancelTime(ridNan) = [];
nStopStop(ridNan) = [];
% if nTrialCriteria == 10
cancelTimeWeightedMean = sum(weightCancelTime .* (nStopStop / sum(nStopStop)))
% elseif nTrialCriteria == 20
%     over20Ind = over20Ind(~nanCancel);
%     nStopStopTotal = sum(nStopStop(over20Ind));
%     cancelTimeWeightedMean = sum(cancelTime(over20Ind) .* (nStopStop(over20Ind) / nStopStopTotal))
% end

%% Mean of cancel times
cancelTimeMean = nanmean(cell2mat(cancelTime))



%% 40 ms ratio around SSRT


% A neuron "cancels" if one of its condition's mean cancel time is within latestTime
goTargSlowSpikeMean = cell(size(cancelData, 1), 1);
stopStopSpikeMean = cell(size(cancelData, 1), 1);
for j = 1 : size(cancelData, 1)
    goTargSlowSpikeMean{j} = cellfun(@mean, cancelData.goTargSlowSpike{j});
    stopStopSpikeMean{j} = cellfun(@mean, cancelData.stopStopSpike{j});
end

ratioAroundSsrt = cellfun(@(in1,in2) in1 ./ in2, goTargSlowSpikeMean, stopStopSpikeMean, 'uni',0);

%% How many (among all units X SSDs were go vs stop different during 40 ms peri-SSRT?
% alphaVal = .05;
%
% peri40msInd = cellfun(@(x) x < alphaVal, cancelData.pValue40msStopStop, 'uni', false);
% notPeri40ms = cellfun(@(x) x >= alphaVal, cancelData.pValue40msStopStop, 'uni', false);
%
% pPeri40ms = sum(cell2mat(peri40msInd)) / (sum(cell2mat(peri40msInd)) + sum(cell2mat(notPeri40ms)));

%% How many ______ had at least one condition with ratio > 1?
ratioLogical = cellfun(@(in1) in1 >=1, ratioAroundSsrt, 'uni', false);


% All units with a condition > 1
nAboveOne = cellfun(@(in1) sum(in1) >=1, ratioLogical);
sum(nAboveOne) / length(nAboveOne)

sum(cell2mat(ratioLogical)) / length(cell2mat(ratioLogical))

%% How many ______ had at least one p value less than .05?
pValueLogical = cellfun(@(in1) in1 <= .05, cancelData.pValue40msStopStop, 'uni', false);

% All units with a condition < .05
disp('How many units at least one below .05')
nBelowCrit = cellfun(@(in1) sum(in1) >= 1, pValueLogical);
sum(nBelowCrit) / length(nBelowCrit)

% All conditions < .05
disp('How many conditions below .05')
sum(cell2mat(pValueLogical)) / length(cell2mat(pValueLogical))

% Mean ratio for those conditions with and without p value less than .05
disp('Stats for those below .05')
belowCritRatio = cell2mat(cellfun(@(in1, in2) in1(in2), ratioAroundSsrt, pValueLogical, 'uni', false));
mean(belowCritRatio)
[h,p,ci,stats] = ttest(belowCritRatio-1)


disp('Stats for those above .05')
aboveCritRatio = cell2mat(cellfun(@(in1, in2) in1(~in2), ratioAroundSsrt, pValueLogical, 'uni', false));
mean(aboveCritRatio)
[h,p,ci,stats] = ttest(aboveCritRatio-1)


