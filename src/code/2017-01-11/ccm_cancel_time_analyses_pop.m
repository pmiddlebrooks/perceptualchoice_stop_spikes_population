function data = ccm_cancel_time_analyses_pop(subject, projectRoot, projectDate, options)

if nargin < 4
    options.multiUnit   = false;
    options.neuronCategory   = 'presacc';
    options.plotFlag   = true;
    options.printPlot   = true;
    options.ssrt = 'intWeightPerSession';
end


dataPath = fullfile(projectRoot,'data',projectDate,subject);

% Plotting variables:
nCol = 2;
nRow = 2;
easyPlot = 1;
hardPlot = 2;

markSizeInd = 30;
markColorInd = [0 0 0];
markSizePop = 120;
markColorPop = [1 0 0];
ylim2Std = [-200 500];

for iCat = 1 : length(options.neuronCategory)
    
    % Build a new table of the relevant neurons, and a list of the session/Unit
    
    % Load the neuron unit list for that category
    fileName = fullfile(dataPath, ['ccm_',options.neuronCategory{iCat},'_neurons']);
    if options.multiUnit
        fileName = [fileName, '_multiUnit'];
    end    
    load(fileName)
    
    % load the population of cancel time anlysis
    canFileName = fullfile(dataPath, 'go_vs_canceled', options.ssrt, ['ccm_canceled_vs_go_neuronTypes']);
    if options.multiUnit
        canFileName = [canFileName, '_multiUnit'];
    end    
    load(canFileName)
    
%     deleteList = ccm_exclude_sessions(subject);
%     deleteInd = ismember(cancelTypes.sessionID, deleteList);
%     cancelTypes(deleteInd, :) = [];
    
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
    
    stopStopEasySsrt = cell(length(ssdArray), 1);
    stopStopHardSsrt = cell(length(ssdArray), 1);
    
    
    if options.plotFlag
        figureHandle = 69;
        xlim2Std = [50 max(ssdArray)+50];
        
        [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(nRow, nCol, figureHandle);
        clf
        
        ax(1, easyPlot) = axes('units', 'centimeters', 'position', [xAxesPosition(1, easyPlot) yAxesPosition(1, easyPlot) axisWidth axisHeight]);
        hold(ax(1, easyPlot), 'on')
        title('2 std neural separation: EASY')
        set(ax(1, easyPlot), 'ylim', ylim2Std, 'xlim', xlim2Std, 'xticklabel', [])
        cla
        
        ax(1, hardPlot) = axes('units', 'centimeters', 'position', [xAxesPosition(1, hardPlot) yAxesPosition(1, hardPlot) axisWidth axisHeight]);
        hold(ax(1, hardPlot), 'on')
        title('2 std neural separation: HARD')
        set(ax(1, hardPlot), 'ylim', ylim2Std, 'xlim', xlim2Std, 'xticklabel', [])
        cla
        
        
        ax(2, easyPlot) = axes('units', 'centimeters', 'position', [xAxesPosition(2, easyPlot) yAxesPosition(2, easyPlot) axisWidth axisHeight]);
        hold(ax(2, easyPlot), 'on')
        title('2std - SSRT neural separation: EASY')
        set(ax(2, easyPlot), 'ylim', ylim2Std, 'xlim', xlim2Std)
        cla
        plot(ax(2, easyPlot), xlim, [0 0], '-k')
        
        ax(2, hardPlot) = axes('units', 'centimeters', 'position', [xAxesPosition(2, hardPlot) yAxesPosition(2, hardPlot) axisWidth axisHeight]);
        hold(ax(2, hardPlot), 'on')
        title('2std - SSRT neural separation: HARD')
        set(ax(2, hardPlot), 'ylim', ylim2Std, 'xlim', xlim2Std)
        cla
        plot(ax(2, hardPlot), xlim, [0 0], '-k')
    end
    
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
            stopStopEasySsrt{i} = cellfun(@(x,y,z) x(y & z), cancelData.stopStopSsrt, easyInd, ssdInd, 'uni', false);
            iUnitEasyInd = cellfun(@(x) ~isempty(x), cancelTimeEasy2Std{i});
            cancelTimeEasySsrt{i} = cell2mat(cancelTimeEasy2Std{i}(iUnitEasyInd)) - cell2mat(stopStopEasySsrt{i}(iUnitEasyInd));
            
            if options.plotFlag
                % 2Std Per neuron per SSD plots
                scatter(ax(1, easyPlot), cell2mat(ssdEasy{i}), cell2mat(cancelTimeEasy2Std{i}), markSizeInd, markColorInd)
                
                % cancel-SSRT Per neuron per SSD plots
                scatter(ax(2, easyPlot), cell2mat(ssdEasy{i}(iUnitEasyInd)), cancelTimeEasySsrt{i}, markSizeInd, markColorInd)
            end
        end
        
        
        if sum(cell2mat(ssdHard{i}))
            cancelTimeHard2Std{i} = cellfun(@(x,y,z,k) x(y & z) - k(y & z), cancelData.cancelTime2Std, hardInd, ssdInd, cancelData.stopStopSsd, 'uni', false);
            stopStopHardSsrt{i} = cellfun(@(x,y,z) x(y & z), cancelData.stopStopSsrt, hardInd, ssdInd, 'uni', false);
            iUnitHardInd = cellfun(@(x) ~isempty(x), cancelTimeHard2Std{i});
            cancelTimeHardSsrt{i} = cell2mat(cancelTimeHard2Std{i}(iUnitHardInd)) - cell2mat(stopStopHardSsrt{i}(iUnitHardInd));
                       
            if options.plotFlag
                % 2Std Per neuron per SSD plots
                scatter(ax(1, hardPlot), cell2mat(ssdHard{i}), cell2mat(cancelTimeHard2Std{i}), markSizeInd, markColorInd)
                
                % cancel-SSRT Per neuron per SSD plots
                scatter(ax(2, hardPlot), cell2mat(ssdHard{i}(iUnitHardInd)), cancelTimeHardSsrt{i}, markSizeInd, markColorInd)
            end
            
        end
        
        
    end
    
    
  
    
    if options.plotFlag
        for i = 1 :length(ssdArray)
            
            if sum(cell2mat(ssdEasy{i}))
                % Per SSD averages
                scatter(ax(1, easyPlot), mean(cell2mat(ssdEasy{i})), nanmean(cell2mat(cancelTimeEasy2Std{i})), markSizePop, 'markerFaceColor', markColorPop, 'markerEdgeColor', 'k')
                
                % cancel-SSRT Per neuron per SSD plots
                scatter(ax(2, easyPlot), mean(cell2mat(ssdEasy{i})), nanmean(cancelTimeEasySsrt{i}), markSizePop, 'markerFaceColor', markColorPop, 'markerEdgeColor', 'k')
            end
            if sum(cell2mat(ssdHard{i}))
                
                % Per SSD averages
                scatter(ax(1, hardPlot), mean(cell2mat(ssdHard{i})), nanmean(cell2mat(cancelTimeHard2Std{i})), markSizePop, 'markerFaceColor', markColorPop, 'markerEdgeColor', 'k')
                
                % cancel-SSRT Per neuron per SSD plots
                scatter(ax(2, hardPlot), mean(cell2mat(ssdHard{i})), nanmean(cancelTimeHardSsrt{i}), markSizePop, 'markerFaceColor', markColorPop, 'markerEdgeColor', 'k')
            end
        end
        
        
        h=axes('Position', [0 0 1 1], 'Visible', 'Off');
        titleString = sprintf('\n%s\tCANCEL TIMES', subject);
        text(0.5,1, titleString, 'HorizontalAlignment','Center', 'VerticalAlignment','Top')
        
        if options.printPlot
            sFileName = fullfile(projectRoot,'results',projectDate,subject,['Cancel_times_population_',options.neuronCategory{iCat}]);
    if options.multiUnit
        sFileName = [sFileName, '_multiUnit'];
    end    
            print(figureHandle, sFileName, '-dpdf', '-r300')
        end
        
    end
end

data.ssdEasy            = ssdEasy;
data.ssdHard            = ssdHard;
data.cancelTimeEasy2Std = cancelTimeEasy2Std;
data.cancelTimeHard2Std = cancelTimeHard2Std;
data.cancelTimeEasySsrt = cancelTimeEasySsrt;
data.cancelTimeHardSsrt = cancelTimeHardSsrt;
data.cancelData         = cancelData; 
