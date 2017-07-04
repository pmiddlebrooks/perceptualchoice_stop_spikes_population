
%%
cancelData.cancelTime = cancelData.cancelTime2Std - cancelData.stopStopSsd - cancelData.ssrt;

test = cellfun(@(x,y,z) x - y - z, cancelData.cancelTime2Std, cancelData.stopStopSsd, num2cell(cancelData.ssrt), 'uni', false)


%%


% load a list of neurons sessions and units
categoryList = {'presacc','presaccNoVis','presaccRamp','visPresacc'};
categoryList = {'presacc'};

% Build a new table of the relevant neurons, and a list of the session/Unit

load(fullfile(dataPath, ['ccm_',categoryList{1},'_neurons']))

% load the population of cancel time anlysis
load(fullfile(dataPath, ['ccm_canceled_vs_go_neuronTypes']))

deleteList = ccm_exclude_sessions(subject);
deleteInd = ismember(cancelTypes.sessionID, deleteList);
cancelTypes(deleteInd, :) = [];

cancelData = table();
for i = 1 : size(neurons, 1)
    % find the indices in cancelTypes that correspond to this unit
    iInd = strcmp(neurons.sessionID(i), cancelTypes.sessionID) & strcmp(neurons.unit(i), cancelTypes.unit);
    cancelData = [cancelData; cancelTypes(iInd,:)];
    
end





ssdArray = unique(cell2mat(cancelData.stopStopSsd));
ssdEasy = cell(length(ssdArray), 1);
ssdHard = cell(length(ssdArray), 1);

cancelTimeEasy2Std = cell(length(ssdArray), 1);
cancelTimeHard2Std = cell(length(ssdArray), 1);

cancelTimeEasy6Std = cell(length(ssdArray), 1);
cancelTimeHard6Std = cell(length(ssdArray), 1);
ssdTimeHard6Std = cell(length(ssdArray), 1);

cancelTimeEasySsrt = cell(length(ssdArray), 1);
cancelTimeHardSsrt = cell(length(ssdArray), 1);
ssdTimeHardSsrt = cell(length(ssdArray), 1);


nCol = 2;
nRow = 2;
easyPlot = 1;
hardPlot = 2;

figureHandle = 68;
[axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(nRow, nCol, figureHandle);
clf

ax(1, easyPlot) = axes('units', 'centimeters', 'position', [xAxesPosition(1, easyPlot) yAxesPosition(1, easyPlot) axisWidth axisHeight]);
hold(ax(1, easyPlot), 'on')
cla
ax(1, hardPlot) = axes('units', 'centimeters', 'position', [xAxesPosition(1, hardPlot) yAxesPosition(1, hardPlot) axisWidth axisHeight]);
hold(ax(1, hardPlot), 'on')
cla


easyInd = cellfun(@(x) x == 1, cancelData.stopStopCond, 'uni', false);
hardInd = cellfun(@(x) x == 2, cancelData.stopStopCond, 'uni', false);

for i = 1 :length(ssdArray)
    ssdInd =  cellfun(@(x) x == ssdArray(i), cancelData.stopStopSsd, 'uni', false);
    
    iEasyCoh = cellfun(@(x,y,z) x(y & z), cancelData.stopStopCond, easyInd, ssdInd, 'uni', false);
    ssdEasy{i} = cellfun(@(x,y,z) x(y & z), cancelData.stopStopSsd, easyInd, ssdInd, 'uni', false);
    
    iHardCoh = cellfun(@(x,y,z) x(y & z), cancelData.stopStopCond, hardInd, ssdInd, 'uni', false);
    ssdHard{i} = cellfun(@(x,y,z) x(y & z), cancelData.stopStopSsd, hardInd, ssdInd, 'uni', false);
    
    fprintf('Easy units: %d\n', sum(~cellfun(@isempty, ssdEasy{i})))
    fprintf('Hard units: %d\n', sum(~cellfun(@isempty, ssdHard{i})))
    if sum(cell2mat(ssdEasy{i}))
        cancelTimeEasy2Std{i} = cellfun(@(x,y,z,k) x(y & z) - k(y & z), cancelData.cancelTime2Std, easyInd, ssdInd, cancelData.stopStopSsd, 'uni', false);
        iUnitEasyInd = cellfun(@(x) ~isempty(x), cancelTimeEasy2Std{i});
        cancelTimeEasySsrt{i} = cell2mat(cancelTimeEasy2Std{i}(iUnitEasyInd)) - cancelData.ssrt(iUnitEasyInd);
        
        scatter(ax(1, easyPlot), cell2mat(ssdEasy{i}), cell2mat(cancelTimeEasy2Std{i}))
        
    end
    if sum(cell2mat(ssdHard{i}))
        cancelTimeHard2Std{i} = cellfun(@(x,y,z,k) x(y & z) - k(y & z), cancelData.cancelTime2Std, hardInd, ssdInd, cancelData.stopStopSsd, 'uni', false);
        iUnitHardInd = cellfun(@(x) ~isempty(x), cancelTimeHard2Std{i});
        cancelTimeHardSsrt{i} = cell2mat(cancelTimeHard2Std{i}(iUnitHardInd)) - cancelData.ssrt(iUnitHardInd);
        
        scatter(ax(1, hardPlot), cell2mat(ssdHard{i}), cell2mat(cancelTimeHard2Std{i}))
    end
    
    
    
    %    hard6Std = cellfun(@(x) ~isnan(x), cancelData.cancelTime6Std, 'uni', false);
    %    cancelTimeEasy6Std{i} =
    
    
    
end


%%


classicEpoch = 'presacc';

% Load Classic neuron classifications
c = load(fullfile(dataPath, ['ccm_',classicEpoch,'_neurons', addMulti]));
classic = c.neurons;

rampEpoch = 'presaccRamp';

% Load Classic neuron classifications
r = load(fullfile(dataPath, ['ccm_',rampEpoch,'_neurons', addMulti]));
ramp = r.neurons;


[~, iCla, iRmp] = setxor(classic, ramp);
classicRamp = intersect(classic, ramp);
classicNoRamp = classic(iCla, :);
rampNoClassic = ramp(iRmp, :);

%%
neuronTypes = classicRamp;
neurons = table();
for i = 61 : size(neuronTypes, 1)
    unitInfo = table();
    unitInfo.sessionID  = neuronTypes.sessionID(i);
    unitInfo.unit       = neuronTypes.unit(i);
    unitInfo.hemisphere  = neuronTypes.hemisphere(i);
    unitInfo.rf         = neuronTypes.rf(i);
    
    fprintf('%02d of %d\t%s\t%s\n',i,size(neuronTypes, 1), neuronTypes.sessionID{i},neuronTypes.unit{i})
    fprintf('Hem: %s\tRF: %s\n',neuronTypes.hemisphere{i},neuronTypes.rf{i})
    
    opt.unitArray = unitInfo.unit;
    opt.hemisphere = neuronTypes.hemisphere{i};
    opt.multiUnit = multiUnit;
    
    pdfName = [neuronTypes.sessionID{i},'_ccm_',neuronTypes.unit{i},'_neuron_collapse.pdf'];
    if exist(fullfile(local_figure_path,subject,'sessionCollapseChoice',pdfName))
        open(fullfile(local_figure_path,subject,'sessionCollapseChoice',pdfName))
    else
        iData = ccm_session_data(subject, neuronTypes.sessionID{i}, opt);
    end
    
    
    prompt = 'add to list?';
    addToList = input(prompt);
    if addToList
        
        neurons = [neurons; unitInfo];
        
    end
    clear iData
end

%%
subject = 'broca';

opt                 = ccm_options;
opt.collapseSignal  = true;
opt.collapseTarget  = true;
opt.doStops         = false;
opt.plotFlag        = false;
opt.printPlot       = false;
opt.figureHandle   	= 106;
opt.multiUnit       = multiUnit;
% opt.multiUnit       = false;

iUnit               = {'bp228n02', 'spikeUnit1'};
iData               = ccm_session_data(subject, iUnit, opt);
%%
iCat              	= ccm_classify_neuron(iData);


%%
psWin = [-149:-99];
psRate = nansum(Data.responseOnset.colorCoh(sigInd).goTarg.raster(:,postsaccAlign + psWin), 2)  .* 1000 ./ length(psWin)

%%

categoryList = {'presaccNoCancel'};
opt = ccm_population_neuron_plot;

opt.doStops = true;
opt.easyOnly = false;
opt.multiUnit = multiUnit;
opt.normalize = true;
opt.categoryName = categoryList{1};
ccm_population_neuron_plot(subject,projectRoot,projectDate,opt)

%%

load(fullfile(dataPath, ['ccm_neuronTypes', addMulti]));

modulated = neuronTypes.fix | neuronTypes.vis | neuronTypes.checker | neuronTypes.presacc | neuronTypes.postsacc;

%%
fileName = fullfile(dataPath, 'go_vs_canceled', ssrtUse, ['ccm_canceled_vs_go_neuronTypes', addMulti]);
load(fileName);
%%
[C,ia,ic] = unique(cancelTypes.sessionID);

cancelTypes = cancelTypes(ia,:);

cancelTypes(cellfun(@isempty, cancelTypes.stopStopSsrt),:) = [];

ssrt = cell2mat(cellfun(@(x) x(1), cancelTypes.stopStopSsrt, 'uni', false))

deleteSessionSsrt = cancelTypes.sessionID(ssrt < 20)
%%

stopStopInd = ismember(cell2table(Data.targOn.easyIn.stopStop.unitStop,'VariableNames',{'sessionID' 'unit'}), neurons(:,1:2))


%%
clf
figure(1)
hold all
i = 1
goInd = find(goTargInd)
plot(Data.checkerOn.easyIn.goTarg.sdf(goInd(i),:))
plot(Data.checkerOn.easyOut.goTarg.sdf(goInd(i),:))

%%

% easyInSpike = mean(Data.checkerOn.easyIn.goTarg.sdf(goTargInd,epochBegin:epochEnd), 2);
% hardInSpike = mean(Data.checkerOn.hardIn.goTarg.sdf(goTargInd,epochBegin:epochEnd), 2);
easyInSpike = sum(Data.checkerOn.easyIn.goTarg.raster(goTargInd,epochBegin:epochEnd), 2) ./ epochDuration;
hardInSpike = sum(Data.checkerOn.hardIn.goTarg.raster(goTargInd,epochBegin:epochEnd), 2) ./ epochDuration;
easyInRT = Data.checkerOn.easyIn.goTarg.rt(goTargInd);
hardInRT = Data.checkerOn.hardIn.goTarg.rt(goTargInd);

rtDiff = easyInRT - hardInRT;
spikeDiff = easyInSpike - hardInSpike;
figure(4)
clf
plot(rtDiff, spikeDiff, '.k')

%%
ssdAll = Data.(opt.epochArray{e}).(conditionArray{cInd}).stopStop.ssd(stopStopInd);
sdfAll = Data.(opt.epochArray{e}).(conditionArray{cInd}).stopStop.sdf(stopStopInd,:);
sdfAllGo = Data.(opt.epochArray{e}).(conditionArray{cInd}).goSlow.sdf(stopStopInd,:);

earlyInd = Data.(opt.epochArray{e}).(conditionArray{cInd}).stopStop.ssd(stopStopInd) < 240;
lateInd = Data.(opt.epochArray{e}).(conditionArray{cInd}).stopStop.ssd(stopStopInd) > 240;

earlySdfMean = nanmean(sdfAll(earlyInd,:), 1);
earlySdfMeanGo = nanmean(sdfAllGo(earlyInd,:), 1);
lateSdfMean = nanmean(sdfAll(lateInd,:), 1);
lateSdfMeanGo = nanmean(sdfAllGo(lateInd,:), 1);

earlySsdMean = nanmean(ssdAll(earlyInd,:));
lateSsdMean = nanmean(ssdAll(lateInd,:));


figure(55)
clf
hold all
plot(earlySdfMean(align+epochRange),'k')
plot(lateSdfMean(align+epochRange),'r')
plot(earlySdfMeanGo(align+epochRange),'g')
plot(lateSdfMeanGo(align+epochRange),'b')

plot(meanGoSlowSdf(align+epochRange),'--k')
plot(meanStopStopSdf(align+epochRange),'--k')

%%
subject = 'joule';

projectDate = '2017-01-11';
projectRoot = '/Volumes/HD-1/Users/paulmiddlebrooks/perceptualchoice_stop_spikes_population';

addpath(genpath(fullfile(projectRoot,'src/code',projectDate)));
dataPath = fullfile(projectRoot,'data',projectDate,subject);


load(fullfile(dataPath, ['ccm_neuronTypes', addMulti]));

sessionList = unique(neuronTypes.sessionID);

nTrial = nan(length(sessionList),1);
nAbort = nan(length(sessionList),1);
for i = 1 : length(sessionList)
    [trialData, SessionData] = load_data(subject, sessionList{i},ccm_min_vars);
    nTrial(i) = length(trialData.rt);
    nAbort(i) = sum(strcmp(trialData.trialOutcome, 'choiceStimulusAbort'));
    
end
sum(nTrial)
mean(nTrial)
sum(nAbort)/sum(nTrial)



%%
subject = 'joule';
% subject = 'broca';

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

%%
dataTable = table();
sessionID = [];
unit = [];
for i = 1 : size(cancelData, 1)
    sessionID = [sessionID; repmat(cancelData.sessionID(i), length(cancelData.nStopStop{i}), 1)];
    unit = [unit; repmat(cancelData.unit(i), length(cancelData.nStopStop{i}), 1)];
end
cancelTime = cellfun(@(x,y,z) x - y - z, cancelData.cancelTime2Std, cancelData.stopStopSsd, cancelData.stopStopSsrt, 'uni', false);


dataTable.sessionID = sessionID;
dataTable.unit = unit;
dataTable.ssrt = cell2mat(cancelData.stopStopSsrt);
dataTable.coherence = cell2mat(cancelData.stopStopCoh);
dataTable.ssd = cell2mat(cancelData.stopStopSsd);
dataTable.nStop = cell2mat(cancelData.nStopStop);
dataTable.pValue40ms = cell2mat(cancelData.pValue40msStopStop);
dataTable.cancelTime = cell2mat(cancelTime);

writetable(dataTable, fullfile(dataPath,'go_vs_canceled',[category,'_cancel_data.csv']))

data = table();
[C,ia,ic] = unique(dataTable.sessionID);
data.sessionID = C;
data.ssrt = dataTable.ssrt(ia);

writetable(data, fullfile(dataPath,'go_vs_canceled',['ssrt_',category,'_cancel_data.csv']))

