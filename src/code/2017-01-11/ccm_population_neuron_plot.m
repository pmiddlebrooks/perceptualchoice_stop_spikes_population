function Data = ccm_population_neuron_plot(subject, projectRoot, projectDate, opt)
%
% function ccm_population_neuron_plot(Data, opt)
%
% Plots avg sdfs from a given population of neurons.
%


%%
% Copied and modified from ccm_session_data_plot
% Set defaults

if nargin < 4
    opt.epochArray      = {'targOn','checkerOn','stopSignalOn','responseOnset'};
    opt.doGos           = true;
    opt.doNoncanceled   = true;
    opt.doCanceled      = true;
    opt.plotError     	= true;
    opt.easyOnly     	= false;
    opt.inOnly          = false;
    opt.normalize         = false;
    opt.excludeSessions = true;
    opt.plotSEM         = false;
    opt.categoryName   	= 'presacc';
    opt.printPlot   	= true;
    opt.dataType        = 'neuron';
    opt.multiUnit     	= false;
    opt.ssrtUse         = 'intWeightPerSession';
    opt.saccadeBaseRatio         = [];
    if nargin == 0, Data = opt; return, end
end

% if ~opt.doStops
%     opt.epochArray      = {'targOn','checkerOn','responseOnset'};
% end

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
printPlot       = true; % opt.printPlot;
% filterData      = opt.filterData;
% stopHz          = opt.stopHz;
% Define axes limits
switch opt.dataType
    case 'neuron'
        if opt.normalize
            yLimMax = 4;
            yLimMin = 1;
        else
            yLimMax = 140;
            yLimMin = 10;
        end
    case {'erp','lfp'}
        yLimMax = .04;
        yLimMin = -.04;
end


goOutcomeArray      = {'goTarg'}; %{'goTarg', 'goDist'};
stopOutcomeArray    = {'stopTarg'}; %{'stopTarg', 'stopDist'};
conditionArray       = {'hardIn', 'hardOut', 'easyIn', 'easyOut'};
conditionArrayInd   = [1 2 3 4];


inStyle = '-';
outStyle = '--';

% goTargStyle = '-';
% goDistStyle = '--';
% stopTargStyle = '-';
% stopDistStyle = '--';
% stopStopStyle = '-';

inOutStyle = {inStyle, outStyle, inStyle, outStyle};  % {'goTarg', 'goDist'};

goEasyColor = [0 .8 0];
goHardColor = [0 .5 0];
if opt.easyOnly, goEasyColor = goHardColor; end
stopEasyColor = [.8 0 0];
stopHardColor = [.5 0 0];
stopStopEasyColor = [.5 .5 .5];
stopStopHardColor = [0 0 0];
% goOutcomeStyle = {goTargStyle, goDistStyle};  % {'goTarg', 'goDist'};
% stopOutcomeStyle = {stopTargStyle, stopDistStyle};  % {'stopTarg', 'stopDist'};

goHardEasyColor = [goHardColor; goHardColor; goEasyColor; goEasyColor];
stopHardEasyColor = [stopHardColor; stopHardColor; stopEasyColor; stopEasyColor];
stopStopColor = [stopStopHardColor; stopStopHardColor;stopStopEasyColor; stopStopEasyColor];
% goAccuracyColor = [goEasyColor; goHardColor];
% stopAccuracyColor = [stopEasyColor; stopHardColor];




if opt.easyOnly
    % conditionArray       = {'easyIn', 'easyOut'};
    conditionArrayInd   = [3 4];
end
if opt.inOnly
    % conditionArray       = {'easyIn', 'hardIn'};
    conditionArrayInd   = [1 3];
end



% ____________________    LOAD DATA    ____________________
dataPath = fullfile(projectRoot,'data',projectDate,subject);

load(fullfile(dataPath, ['ccm_',opt.categoryName,'_neurons',addMulti])) % Load Data struct for that population
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

Data = load(fullfile(dataPath, ['ccm_all_neuron_population',addMulti,addNorm])); % Load Data struct for that population



% load the population of cancel time anlysis
if opt.doCanceled
fileName = fullfile(dataPath, 'go_vs_canceled', opt.ssrtUse, ['ccm_canceled_vs_go_neuronTypes', addMulti]);
elseif opt.doNoncanceled
fileName = fullfile(dataPath, 'go_vs_noncanceled', opt.ssrtUse, ['ccm_noncanceled_vs_go_neuronTypes', addMulti]);
end
load(fileName)









% %  GET STOPPING/CANCELING POPULATION Data
% 
% 
% % Index the relelvant units for the StopStop Data
% unitIndEasy = ismember(cell2table(Data.checkerOn.easyIn.stopStop.unitStop,'VariableNames',{'sessionID' 'unit'}), neurons(:,1:2));
% unitIndHard = ismember(cell2table(Data.checkerOn.hardIn.stopStop.unitStop,'VariableNames',{'sessionID' 'unit'}), neurons(:,1:2));
% % Concatenate all data from the same units
% cancelEasy = cell2table(Data.checkerOn.easyIn.stopStop.unitStop(unitIndEasy,:),'VariableNames',{'sessionID' 'unit'});
% cancelEasy.ssd = Data.checkerOn.easyIn.stopStop.ssd(unitIndEasy);
% cancelHard = cell2table(Data.checkerOn.hardIn.stopStop.unitStop(unitIndHard,:),'VariableNames',{'sessionID' 'unit'});
% cancelHard.ssd = Data.checkerOn.hardIn.stopStop.ssd(unitIndHard);
% 
% 
% % Build a table of ssrt and cancel time
% cancelEasy.ssrt = nan(length(cancelEasy.ssd), 1);
% cancelEasy.cancelTime = nan(length(cancelEasy.ssd), 1);
% cancelHard.ssrt = nan(length(cancelHard.ssd), 1);
% cancelHard.cancelTime = nan(length(cancelHard.ssd), 1);
% 
% for i = 1 : length(cancelEasy.ssd)
%     iUnitInd = ismember(cancelTypes.sessionID, cancelEasy.sessionID(i)) & ismember(cancelTypes.unit, cancelEasy.unit(i));
%     iSsdInd = cancelTypes.stopStopSsd{iUnitInd} == cancelEasy.ssd(i);
%     iCondInd = cancelTypes.stopStopCond{iUnitInd} == 1;
%     cancelEasy.ssrt(i) = cancelTypes.stopStopSsrt{iUnitInd}(iSsdInd & iCondInd);
%     cancelEasy.cancelTime(i) = cancelTypes.cancelTime2Std{iUnitInd}(iSsdInd & iCondInd) - cancelEasy.ssrt(i) - cancelEasy.ssd(i);
% end
% for i = 1 : length(cancelHard.ssd)
%     iUnitInd = ismember(cancelTypes.sessionID, cancelHard.sessionID(i)) & ismember(cancelTypes.unit, cancelHard.unit(i));
%     iSsdInd = cancelTypes.stopStopSsd{iUnitInd} == cancelHard.ssd(i);
%     iCondInd = cancelTypes.stopStopCond{iUnitInd} == 2;
%     cancelHard.ssrt(i) = cancelTypes.stopStopSsrt{iUnitInd}(iSsdInd & iCondInd);
%     cancelHard.cancelTime(i) = cancelTypes.cancelTime2Std{iUnitInd}(iSsdInd & iCondInd) - cancelHard.ssrt(i) - cancelHard.ssd(i);
% end











%   ____________________ SET UP PLOT  ____________________

lineWidth = 4;  % for all conditions right now

nCol = length(opt.epochArray);
nRow = 1;

figureHandle = 846;


if printPlot
    [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(nRow, nCol, figureHandle);
else
    [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = screen_figure(nRow, nCol, figureHandle);
end
clf


% _______  Set up axes  ___________
% axes names
axGo = 1;
axStop = 1;
axStop = 1;
for e = 1 : length(opt.epochArray)
    epochRange = ccm_epoch_range(opt.epochArray{e}, 'plot');
    
    
    % Set up plots
    % Go trials
    ax(axGo, e) = axes('units', 'centimeters', 'position', [xAxesPosition(axGo, e) yAxesPosition(axGo, e) axisWidth axisHeight]);
    hold(ax(axGo, e), 'on')
    %             set(ax(axGo, e), 'ylim', [yLimMin yLimMax], 'xlim', [epochRange(1) epochRange(end)])
    set(ax(axGo, e), 'ylim', [yLimMin yLimMax], 'xlim', [1 epochRange(end) - epochRange(1)])
    cla
    hold(ax(axGo, e), 'on')
    plot(ax(axGo, e), [-epochRange(1) -epochRange(1)], [yLimMin yLimMax * .9], '-k', 'linewidth', 2)
    title(opt.epochArray{e})
    set(ax(axGo, e), 'xticklabel', get(gca, 'xtick')+epochRange(1)-1)
    if e > 1, set(ax(axGo, e), 'yticklabel', []), end
    
end







% Loop through epochs
for e = 1 : length(opt.epochArray)
    epochRange = ccm_epoch_range(opt.epochArray{e}, 'plot');
    
    % Loop colorherence (easyIn, easyOut, hardIn, hardOut)
    for c = 1 : length(conditionArrayInd)
        cInd = conditionArrayInd(c);
        
        
        
        
        
            
            
            %   _______ STOP TARG TRIALS  _______
            if opt.doNoncanceled
                % Use the neurons table to Index the relelvant units from the Data struct for StopStop
                stopTargInd = ismember(cell2table(Data.(opt.epochArray{e}).(conditionArray{cInd}).stopTarg.unitStop,'VariableNames',{'sessionID' 'unit'}), neurons(:,1:2));
                meanStopTargSdf = nanmean(Data.(opt.epochArray{e}).(conditionArray{cInd}).stopTarg.sdf(stopTargInd,:), 1);
                if ~strcmp(opt.epochArray{e}, 'stopSignalOn')
                    meanGoFastSdf = nanmean(Data.(opt.epochArray{e}).(conditionArray{cInd}).goFast.sdf(stopTargInd,:), 1);
                end
                
                
                % Get SSRT for the session
                if strcmp(opt.epochArray{e}, 'checkerOn')
                    meanStopTargSsd = nanmean(Data.(opt.epochArray{e}).(conditionArray{cInd}).stopTarg.ssd(stopTargInd), 1);
                    
                    stopTargSsrt = [];
                    
                    for i = find(stopTargInd)'
                        % What's the cancelTypes index for this unit?
                        iCancelTypesInd = strcmp(Data.(opt.epochArray{e}).(conditionArray{cInd}).stopTarg.unitStop(i,1), noncancelTypes.sessionID) & strcmp(Data.(opt.epochArray{e}).(conditionArray{cInd}).stopTarg.unitStop(i,2), noncancelTypes.unit);
                        
                        % What's the SSD/SSRT/Cancel time index (for this unit in the noncancelTypes table)
                        switch conditionArray{cInd}
                            case {'easyIn', 'easyOut'}
                                iCond = 1;
                            case {'hardIn', 'hardOut'}
                                iCond = 2;
                        end
                        iSsdInd = noncancelTypes.stopTargSsd{iCancelTypesInd} == Data.(opt.epochArray{e}).(conditionArray{cInd}).stopTarg.ssd(i) & ...
                            noncancelTypes.stopTargCond{iCancelTypesInd} == iCond;
                        
                        stopTargSsrt = [stopTargSsrt;  noncancelTypes.stopTargSsrt{iCancelTypesInd}(iSsdInd)];
                        
                    end
                end
                
                
                              
                
                align = Data.(opt.epochArray{e}).alignStop;
                if strcmp(opt.epochArray{e}, 'checkerOn')
                    meanStopTargSsrt = nanmean(stopTargSsrt, 1);
                    
%                     plot(ax(axStop, 2), [meanStopTargSsd + -epochRange(1), meanStopTargSsd + -epochRange(1)] , [yLimMin yLimMax * .7], '-', 'color', stopHardEasyColor(cInd,:), 'linewidth', 2)
%                     plot(ax(axStop, 2), [meanStopTargSsrt + meanStopTargSsd + -epochRange(1), meanStopTargSsrt + meanStopTargSsd + -epochRange(1)] , [yLimMin yLimMax * .7], '--', 'color', stopHardEasyColor(cInd,:), 'linewidth', 2)
                end
                
                plot(ax(axStop, e), meanStopTargSdf(align+epochRange(1):end), 'LineWidth',lineWidth,'LineStyle',inOutStyle{cInd},'color',stopHardEasyColor(cInd,:))
                if ~strcmp(opt.epochArray{e}, 'stopSignalOn')
                    plot(ax(axStop, e), meanGoFastSdf(align+epochRange(1):end), 'LineWidth',lineWidth,'LineStyle',inOutStyle{cInd},'color',goHardEasyColor(cInd,:))
                end
                if e == 1
                    fprintf('StopTarg %s conditions: %d\n', conditionArray{cInd}, sum(stopTargInd));
                end
                
                %         if ~strcmp(opt.epochArray{e}, 'stopSignalOn')
                %             meanSDF = nanmean(Data.(conditionArray{cInd}).goFast.(opt.epochArray{e}).sdf, 1);
                %             align = Data.(conditionArray{cInd}).goFast.(opt.epochArray{e}).align;
                %             plot(ax(axStop, e), meanSDF(align+epochRange(1):end), 'LineWidth',lineWidth,'LineStyle',inOutStyle{cInd},'color',goHardEasyColor(cInd,:))
                %         end
            end
            
            
            
            
            
            
%             %   _______ STOP STOP TRIALS  _______
%             
%             if opt.doCanceled
%                 
%                 % Use the neurons table to Index the relelvant units from the Data struct for StopStop
%                 stopStopInd = ismember(cell2table(Data.(opt.epochArray{e}).(conditionArray{cInd}).stopStop.unitStop,'VariableNames',{'sessionID' 'unit'}), neurons(:,1:2));
%                 meanStopStopSdf = nanmean(Data.(opt.epochArray{e}).(conditionArray{cInd}).stopStop.sdf(stopStopInd,:), 1);
%                 if ~strcmp(opt.epochArray{e}, 'stopSignalOn')
%                     meanGoSlowSdf = nanmean(Data.(opt.epochArray{e}).(conditionArray{cInd}).goSlow.sdf(stopStopInd,:), 1);
%                 end
%                 
%                 % Get SSRT for the session, and cancel time for the unit/condition
%                 if strcmp(opt.epochArray{e}, 'checkerOn')
%                     meanStopStopSsd = nanmean(Data.(opt.epochArray{e}).(conditionArray{cInd}).stopStop.ssd(stopStopInd), 1);
%                     
%                     stopStopSsrt = [];
%                     stopStopCancel = [];
%                     
%                     for i = find(stopStopInd)'
%                         % What's the cancelTypes index for this unit?
%                         iCancelTypesInd = strcmp(Data.(opt.epochArray{e}).(conditionArray{cInd}).stopStop.unitStop(i,1), cancelTypes.sessionID) & strcmp(Data.(opt.epochArray{e}).(conditionArray{cInd}).stopStop.unitStop(i,2), cancelTypes.unit);
%                         
%                         % What's the SSD/SSRT/Cancel time index (for this unit in the cancelTypes table)
%                         switch conditionArray{cInd}
%                             case {'easyIn', 'easyOut'}
%                                 iCond = 1;
%                             case {'hardIn', 'hardOut'}
%                                 iCond = 2;
%                         end
%                         iSsdInd = cancelTypes.stopStopSsd{iCancelTypesInd} == Data.(opt.epochArray{e}).(conditionArray{cInd}).stopStop.ssd(i) & ...
%                             cancelTypes.stopStopCond{iCancelTypesInd} == iCond;
%                         
%                         stopStopSsrt = [stopStopSsrt;  cancelTypes.stopStopSsrt{iCancelTypesInd}(iSsdInd)];
%                         stopStopCancel = [stopStopCancel; cancelTypes.cancelTime2Std{iCancelTypesInd}(iSsdInd)];
%                         
%                         
%                         
%                         
%                         
%                         
%                         %                         figure(77);
%                         %                         clf
%                         %                         hold 'all'
%                         %                         iGoSlowSdf = Data.(opt.epochArray{e}).(conditionArray{cInd}).goSlow.sdf(i,:);
%                         %                         iStopStopSdf = Data.(opt.epochArray{e}).(conditionArray{cInd}).stopStop.sdf(i,:);
%                         %                         iSsd = Data.(opt.epochArray{e}).(conditionArray{cInd}).stopStop.ssd(i);
%                         %                         iSsrt = cancelTypes.stopStopSsrt{iCancelTypesInd}(iSsdInd);
%                         %                         iCancelTime = cancelTypes.cancelTime2Std{iCancelTypesInd}(iSsdInd);
%                         %
%                         %                     plot(iGoSlowSdf(align+epochRange(1):end), 'LineWidth',lineWidth,'LineStyle',inOutStyle{cInd},'color',goHardEasyColor(cInd,:))
%                         %                     plot(iStopStopSdf(align+epochRange(1):end), 'LineWidth',lineWidth,'LineStyle',inOutStyle{cInd},'color',stopStopColor(cInd,:))
%                         %                     plot([iSsd + -epochRange(1), iSsd + -epochRange(1)] , [yLimMin yLimMax * .7], '-', 'color', stopStopColor(cInd,:), 'linewidth', 2)
%                         %                     plot([iSsrt + iSsd + -epochRange(1), iSsrt + iSsd + -epochRange(1)] , [yLimMin yLimMax * .7], '--', 'color', stopStopColor(cInd,:), 'linewidth', 2)
%                         %                     plot([iCancelTime + -epochRange(1), iCancelTime + -epochRange(1)] , [yLimMin yLimMax * .7], ':', 'color', stopStopColor(cInd,:), 'linewidth', 2)
%                         %
%                         %     print(77, fullfile(local_figure_path,subject,'cancel_examples',[opt.categoryName,'_',Data.(opt.epochArray{e}).(conditionArray{cInd}).stopStop.unitStop{i,1},'_',Data.(opt.epochArray{e}).(conditionArray{cInd}).stopStop.unitStop{i,2},'_',conditionArray{cInd},'_',num2str(iSsd)]), '-dpdf', '-r300')
%                         
%                         
%                         
%                         
%                     end
%                     
%                 end
%                 
%                 
%                 align = Data.(opt.epochArray{e}).alignStop;
%                 
%                 if strcmp(opt.epochArray{e}, 'checkerOn')
%                     meanStopStopSsrt = nanmean(stopStopSsrt, 1);
%                     meanStopStopCancel = nanmean(stopStopCancel, 1);
%                     
%                     plot(ax(axStop, 2), [meanStopStopSsd + -epochRange(1), meanStopStopSsd + -epochRange(1)] , [yLimMin yLimMax * .7], '-', 'color', stopStopColor(cInd,:), 'linewidth', 2)
%                     plot(ax(axStop, 2), [meanStopStopSsrt + meanStopStopSsd + -epochRange(1), meanStopStopSsrt + meanStopStopSsd + -epochRange(1)] , [yLimMin yLimMax * .7], '--', 'color', stopStopColor(cInd,:), 'linewidth', 2)
%                     plot(ax(axStop, 2), [meanStopStopCancel + -epochRange(1), meanStopStopCancel + -epochRange(1)] , [yLimMin yLimMax * .7], ':', 'color', stopStopColor(cInd,:), 'linewidth', 2)
%                 end
%                 
%                 %                 if ~strcmp(opt.epochArray{e}, 'responseOnset')
%                 plot(ax(axStop, e), meanStopStopSdf(align+epochRange(1):end), 'LineWidth',lineWidth,'LineStyle',inOutStyle{cInd},'color',stopStopColor(cInd,:))
%                 %                 end
%                 if ~strcmp(opt.epochArray{e}, 'stopSignalOn')
%                     plot(ax(axStop, e), meanGoSlowSdf(align+epochRange(1):end), 'LineWidth',lineWidth,'LineStyle',inOutStyle{cInd},'color',goHardEasyColor(cInd,:))
%                 end
%                 if e == 1
%                     fprintf('StopStop %s conditions: %d\n', conditionArray{cInd}, sum(stopStopInd));
%                 end
%                 
%                 %         if ~strcmp(opt.epochArray{e}, 'stopSignalOn')
%                 %             meanSDF = nanmean(Data.(conditionArray{cInd}).goSlow.(opt.epochArray{e}).sdf, 1);
%                 %             align = Data.(conditionArray{cInd}).goSlow.(opt.epochArray{e}).align;
%                 %             plot(ax(axStop, e), meanSDF(align+epochRange(1):end), 'LineWidth',lineWidth,'LineStyle',inOutStyle{cInd},'color',goInOutColor(cInd,:))
%                 %         end
%                 
%             end
%             
            
        



        
        
            %   _______ STOP STOP TRIALS  _______
            
            if opt.doCanceled
                
                % Use the neurons table to Index the relelvant units from the Data struct for StopStop
                stopStopEarlyInd = ismember(cell2table(Data.(opt.epochArray{e}).(conditionArray{cInd}).stopStop.unitStop,'VariableNames',{'sessionID' 'unit'}), neurons(:,1:2));
                
%                 ssdAll = unique(Data.(opt.epochArray{e}).(conditionArray{cInd}).stopStop.ssd(stopStopInd));
%                 stopStopEarly = ssdAll(ceil(length(ssdAll)/2));
%                 stopStopEarly = ssdAll(ceil(2*length(ssdAll)/3));
%                 stopStopEarlyInd = stopStopInd & Data.(opt.epochArray{e}).(conditionArray{cInd}).stopStop.ssd < stopStopEarly;
                
                meanStopStopSdf = nanmean(Data.(opt.epochArray{e}).(conditionArray{cInd}).stopStop.sdf(stopStopEarlyInd,:), 1);
                if ~strcmp(opt.epochArray{e}, 'stopSignalOn')
                    meanGoSlowSdf = nanmean(Data.(opt.epochArray{e}).(conditionArray{cInd}).goSlow.sdf(stopStopEarlyInd,:), 1);
                end
                
                % Get SSRT for the session, and cancel time for the unit/condition
                if strcmp(opt.epochArray{e}, 'checkerOn')
                    meanStopStopSsd = nanmean(Data.(opt.epochArray{e}).(conditionArray{cInd}).stopStop.ssd(stopStopEarlyInd), 1);
                    
                    stopStopSsrt = [];
                    stopStopCancel = [];
                    
                    for i = find(stopStopEarlyInd)'
                        % What's the cancelTypes index for this unit?
                        iCancelTypesInd = strcmp(Data.(opt.epochArray{e}).(conditionArray{cInd}).stopStop.unitStop(i,1), cancelTypes.sessionID) & strcmp(Data.(opt.epochArray{e}).(conditionArray{cInd}).stopStop.unitStop(i,2), cancelTypes.unit);
                        
                        % What's the SSD/SSRT/Cancel time index (for this unit in the cancelTypes table)
                        switch conditionArray{cInd}
                            case {'easyIn', 'easyOut'}
                                iCond = 1;
                            case {'hardIn', 'hardOut'}
                                iCond = 2;
                        end
                        iSsdInd = cancelTypes.stopStopSsd{iCancelTypesInd} == Data.(opt.epochArray{e}).(conditionArray{cInd}).stopStop.ssd(i) & ...
                            cancelTypes.stopStopCond{iCancelTypesInd} == iCond;
                        
                        stopStopSsrt = [stopStopSsrt;  cancelTypes.stopStopSsrt{iCancelTypesInd}(iSsdInd)];
                        stopStopCancel = [stopStopCancel; cancelTypes.cancelTime2Std{iCancelTypesInd}(iSsdInd)];
                        
                        
                        
                        
                        
                        
                        %                         figure(77);
                        %                         clf
                        %                         hold 'all'
                        %                         iGoSlowSdf = Data.(opt.epochArray{e}).(conditionArray{cInd}).goSlow.sdf(i,:);
                        %                         iStopStopSdf = Data.(opt.epochArray{e}).(conditionArray{cInd}).stopStop.sdf(i,:);
                        %                         iSsd = Data.(opt.epochArray{e}).(conditionArray{cInd}).stopStop.ssd(i);
                        %                         iSsrt = cancelTypes.stopStopSsrt{iCancelTypesInd}(iSsdInd);
                        %                         iCancelTime = cancelTypes.cancelTime2Std{iCancelTypesInd}(iSsdInd);
                        %
                        %                     plot(iGoSlowSdf(align+epochRange(1):end), 'LineWidth',lineWidth,'LineStyle',inOutStyle{cInd},'color',goHardEasyColor(cInd,:))
                        %                     plot(iStopStopSdf(align+epochRange(1):end), 'LineWidth',lineWidth,'LineStyle',inOutStyle{cInd},'color',stopStopColor(cInd,:))
                        %                     plot([iSsd + -epochRange(1), iSsd + -epochRange(1)] , [yLimMin yLimMax * .7], '-', 'color', stopStopColor(cInd,:), 'linewidth', 2)
                        %                     plot([iSsrt + iSsd + -epochRange(1), iSsrt + iSsd + -epochRange(1)] , [yLimMin yLimMax * .7], '--', 'color', stopStopColor(cInd,:), 'linewidth', 2)
                        %                     plot([iCancelTime + -epochRange(1), iCancelTime + -epochRange(1)] , [yLimMin yLimMax * .7], ':', 'color', stopStopColor(cInd,:), 'linewidth', 2)
                        %
                        %     print(77, fullfile(local_figure_path,subject,'cancel_examples',[opt.categoryName,'_',Data.(opt.epochArray{e}).(conditionArray{cInd}).stopStop.unitStop{i,1},'_',Data.(opt.epochArray{e}).(conditionArray{cInd}).stopStop.unitStop{i,2},'_',conditionArray{cInd},'_',num2str(iSsd)]), '-dpdf', '-r300')
                        
                        
                        
                        
                    end
                    
                end
                
                
                align = Data.(opt.epochArray{e}).alignStop;
                
                if strcmp(opt.epochArray{e}, 'checkerOn')
                    meanStopStopSsrt = nanmean(stopStopSsrt, 1);
                    meanStopStopCancel = nanmean(stopStopCancel, 1);
                    
                    plot(ax(axStop, 2), [meanStopStopSsd + -epochRange(1), meanStopStopSsd + -epochRange(1)] , [yLimMin yLimMax * .7], '-', 'color', stopStopColor(cInd,:), 'linewidth', 2)
                    plot(ax(axStop, 2), [meanStopStopSsrt + meanStopStopSsd + -epochRange(1), meanStopStopSsrt + meanStopStopSsd + -epochRange(1)] , [yLimMin yLimMax * .7], '--', 'color', stopStopColor(cInd,:), 'linewidth', 2)
                    plot(ax(axStop, 2), [meanStopStopCancel + -epochRange(1), meanStopStopCancel + -epochRange(1)] , [yLimMin yLimMax * .7], ':', 'color', stopStopColor(cInd,:), 'linewidth', 2)
                end
                
                %                 if ~strcmp(opt.epochArray{e}, 'responseOnset')
                plot(ax(axStop, e), meanStopStopSdf(align+epochRange(1):end), 'LineWidth',lineWidth,'LineStyle',inOutStyle{cInd},'color',stopStopColor(cInd,:))
                %                 end
                if ~strcmp(opt.epochArray{e}, 'stopSignalOn')
                    plot(ax(axStop, e), meanGoSlowSdf(align+epochRange(1):end), 'LineWidth',lineWidth,'LineStyle',inOutStyle{cInd},'color',goHardEasyColor(cInd,:))
                end
                if e == 1
                    fprintf('StopStop %s conditions: %d\n', conditionArray{cInd}, sum(stopStopEarlyInd));
                end
                
                %         if ~strcmp(opt.epochArray{e}, 'stopSignalOn')
                %             meanSDF = nanmean(Data.(conditionArray{cInd}).goSlow.(opt.epochArray{e}).sdf, 1);
                %             align = Data.(conditionArray{cInd}).goSlow.(opt.epochArray{e}).align;
                %             plot(ax(axStop, e), meanSDF(align+epochRange(1):end), 'LineWidth',lineWidth,'LineStyle',inOutStyle{cInd},'color',goInOutColor(cInd,:))
                %         end
                
            end
            

            
            
            
            
            
            
            
            
            
            
            
            
            
        
        if opt.doGos
            %   _______ GO TRIALS  _______
% Index the relelvant units for the Go Data
goTargInd = ismember(cell2table(Data.unitGo,'VariableNames',{'sessionID' 'unit'}), neurons(:,1:2));

if ~strcmp(opt.epochArray{e}, 'stopSignalOn')
                % Go Targ
                meanSDF = nanmean(Data.(opt.epochArray{e}).(conditionArray{cInd}).goTarg.sdf(goTargInd,:), 1);
                %             align = Data.alignGo(e);
                align = Data.(opt.epochArray{e}).alignGo;
                plot(ax(axGo, e), meanSDF(align+epochRange(1):end), 'LineWidth',lineWidth,'LineStyle',inOutStyle{cInd},'color',goHardEasyColor(cInd,:))
                if e == 1
                    fprintf('GoTarg %s conditions: %d\n', conditionArray{cInd}, sum(goTargInd));
                end
                
                %                 meanSDF = nanmean(Data.(conditionArray{cInd}).goDist.(opt.epochArray{e}).sdf, 1);
                %                 align = Data.(conditionArray{cInd}).goDist.(opt.epochArray{e}).align;
                %                 plot(ax(axGo, e), meanSDF(align+epochRange(1):end), 'LineWidth',lineWidth,'LineStyle',inOutStyle{cInd},'color',goInOutColor(cInd,:))
            end
        end
        
        
    end % colorCohArray
    
end % epochs


h=axes('Position', [0 0 1 1], 'Visible', 'Off');
titleString = sprintf('%s\t n = %d', opt.categoryName , size(neurons, 1));
text(0.5,1, titleString, 'HorizontalAlignment','Center', 'VerticalAlignment','Top')



if opt.printPlot
    if opt.doGos
        appendFile = '_Go';
    elseif opt.doCanceled
        appendFile = '_Canceled';
    elseif opt.doNoncanceled
        appendFile = '_Noncanceled';
    end
    
    if ~isempty(opt.saccadeBaseRatio)
        appendBase = '_saccBaseRatio_';
    else
        appendBase = [];
    end
    
    filePath = fullfile(projectRoot,'results',projectDate,subject);
    if opt.easyOnly
        fileName = ['pop_',opt.categoryName,'_easy',addMulti,addNorm,appendBase, appendFile,'.pdf'];
    elseif opt.doGos
        fileName = ['pop_',opt.categoryName,'_go',addMulti,addNorm,appendBase, appendFile,'.pdf'];
    else
        fileName = ['pop_',opt.categoryName,addMulti,addNorm,appendBase, appendFile,'.pdf'];
    end
    print(figureHandle, fullfile(filePath, fileName), '-dpdf', '-r300')
end
return

% ________________________________________________________________________________
% ________________________________________________________________________________
% ________________________________________________________________________________
% ________________________________________________________________________________
% ________________________________________________________________________________



