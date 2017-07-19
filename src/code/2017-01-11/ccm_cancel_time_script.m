%%

subject = 'joule';
subject = 'broca';


multiUnit = true;
deleteUnmodulated = true;
deleteSessions = true;
normalizeData = false;
append = false;
ssrtUse = 'intWeightPerSession';

if multiUnit
    addMulti = '_multiUnit';
else
    addMulti = [];
end
if normalizeData
    addNorm = '_normalized';
else
    addNorm = [];
end



sessionRemove = ccm_exclude_sessions(subject);

saccadeBaseRatio = [];
saccadeBaseRatio = 2;

category = 'presacc_cancel_meanDifference';
% category = 'presacc_ddmRankMeanStim_cancel_meanDifference';
% category = 'presacc_cancel_trialByTrial';
% category = 'presacc_ddmRankMeanStim_cancel_trialByTrial';
% category = 'presacc_cancel_meanSdf';
% category = 'presacc_ddmRankMeanStim_cancel_meanSdf';

nTrialCriteria = 10;
% nTrialCriteria = 20;
% nTrialCriteria = 1020;




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
%     find(iInd)
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
% cancelTime = cellfun(@(x,y,z) x - y - z, cancelData.cancelTime2Std, cancelData.stopStopSsd, cancelData.stopStopSsrt, 'uni', false);

% Grand average cancel time
% ssrtGrandAvg = mean(cell2mat(cancelData.stopStopSsrt));
% cancelTime = cellfun(@(x,y) x - y - ssrtGrandAvg, cancelData.cancelTime2Std, cancelData.stopStopSsd, 'uni', false);
% 

cancelTime = cellfun(@(x,y) x - y - 82, cancelData.cancelTime2Std, cancelData.stopStopSsd, 'uni', false);

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

%% Mean of cancel times
cancelTimeMean = nanmean(cell2mat(cancelTime))



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










%% Cancel time using peri-SSD as baseline

% Make cancelData table more convenient by reassigning it to "c".
c = cancelData;

% Constants
% -------------------

ms2Std = 30;
window = 0 : 999; % define an epoch of analysis
preBase = 20; % ms before SSD to begin baseline epoch
postBase = 40; % ms after SSD to end baseline epoch

% Initialize cancelTime to add to c table
cancelTimeNewBase = cell(size(c, 1));

% Loop through each row (unit) of c table
for i = 1 : size(c, 1)
    nCond = length(c.stopStopSsrt(i));
    iCancelTime = nan(nCond, 1);
    
    % Loop through each valid condition of each unit
    for j = 1 : nCond
        jCoh = c.stopStopCoh{i}(j);
        jSsd = c.stopStopSsd{i}(j);
        
        % Get the differential SDF to analyze
        jSdfDiff = c.goTargSlowSpike{i}(j,c.goTargSlowCheckerAlign{i}(j)+window) - c.stopStopSpike{i}(j,c.stopStopCheckerAlign{i}(j)+window);
        
        % Use the baseline epoch to calculate 2 Std of differential SDF
        j2Std = 2*std(jSdfDiff(jSsd-preBase : jSsd+postBase));
        
                % are there times at which the difference between sdfs is
                % greater than 2 standard deviations of the differential
                % SDF around the SSD?
                std2Ind = jSdfDiff(jSsd:end) > j2Std;

                
                
                
                
                
                
                                % Look for a sequence of ms2Std ms for which the go sdf is 2
                % std greater than the stop sdf.
                % First whether the differential sdf was > 2*Std for the
                % first ms2Std ms
                if sum(std2Ind(1:ms2Std)) == ms2Std
                    cancelTimeNewBase{i}(j) = c.stopStopSsd{i}(j);
                else
                    % If it wasn't, detmerine whether there was a time
                    % after SSD that the differential
                    % sdf was > 2*Std for at least ms2Std ms.
                    riseAbove2Std = find([0; diff(std2Ind)] == 1);
                    sinkBelow2Std = find([0; diff(std2Ind)] == -1);
                    if ~isempty(riseAbove2Std)                        
                        % If there's one more riseAbove2Std than sinkBelow2Std, the last riseAbove2Std
                        % will last until the end of the sdf: Add to
                        % singkBelowStd the end of the epoch
                        if length(riseAbove2Std) > length(sinkBelow2Std)
                            sinkBelow2Std = [sinkBelow2Std; checkerEpochEnd];
                        end
                        
                        % Now riseAbove2Std length should be equal. See if
                        % any of the riseAbove2Std streaks go longer than
                        % 50ms
                        ind = find(sinkBelow2Std - riseAbove2Std >= ms2Std, 1);
                        if ~isempty(ind)
                            cancelTimeNewBase{i}(j) = riseAbove2Std(ind) + stopStopSsd(i);
                        end
                     end
                end

                
                
                
                
        
    end
end

dataTable.sessionID     = sessionID;
dataTable.unit          = unit;
dataTable.ssrt          = cell2mat(c.stopStopSsrt);
dataTable.coherence     = cell2mat(c.stopStopCoh);
dataTable.ssd           = cell2mat(c.stopStopSsd);
dataTable.nStop         = cell2mat(c.nStopStop);
dataTable.pValue40ms    = cell2mat(c.pValue40msStopStop);
dataTable.cancelTime    = cell2mat(cancelTime);


% writetable(dataTable, fullfile(dataPath,'go_vs_canceled',['ssrt_',category,'_cancel_data.csv']))

