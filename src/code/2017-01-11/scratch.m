
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





