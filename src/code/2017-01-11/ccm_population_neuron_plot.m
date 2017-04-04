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
    opt.epochArray      = {'targOn','checkerOn','stopSignalOn','responseOnset','rewardOn'};
    opt.epochArray      = {'targOn','checkerOn','stopSignalOn','responseOnset'};
    opt.doStops         = true;
    opt.plotError     	= true;
    opt.easyOnly     	= false;
    opt.inOnly          = false;
    opt.normalize         = false;
    opt.plotSEM         = false;
    opt.categoryName   	= 'presacc';
    opt.printPlot   	= true;
    opt.dataType        = 'neuron';
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
printPlot       = true; % opt.printPlot;
% filterData      = opt.filterData;
% stopHz          = opt.stopHz;
% Define axes limits
switch opt.dataType
    case 'neuron'
        if opt.normalize
            yLimMax = 1;
            yLimMin = .2;
        else
            yLimMax = 110;
            yLimMin = 20;
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

goInOutColor = [goHardColor; goHardColor; goEasyColor; goEasyColor];
stopInOutColor = [stopHardColor; stopHardColor; stopEasyColor; stopEasyColor];
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
Data = load(fullfile(dataPath, ['ccm_all_neuron_population',addMulti,addNorm])); % Load Data struct for that population

% Index the relelvant units for the Go Data
goTargInd = ismember(cell2table(Data.unitGo,'VariableNames',{'sessionID' 'unit'}), neurons(:,1:2));


% if opt.normalizeData
%
%
%     nUnitGo = length(Data.checkerOn.easyIn.goTarg.rt);
%     for i = 1 : nSession
%         iNormFactor = max(Data.responseOnset.easyIn.goTarg.sdf(i,:));
%
%         for e = 1 : length(opt.epochArray)
%             for c = 1 : length(conditionArray)
%                 % Go Data
%                 for g = 1 : length(goOutcomeArray)
%                     if ~strcmp(opt.epochArray{e}, 'stopSignalOn')
%
%                         Data.(opt.epochArray{e}).(conditionArray{c}).(goOutcomeArray{g}).sdf(i,:) = ...
%                             Data.(opt.epochArray{e}).(conditionArray{c}).(goOutcomeArray{g}).sdf(i,:) / iNormFactor;
%                     end
%                 end
%
%                  % Stop Data
%                for s = 1 : length(stopOutcomeArray)
%                    i
%                    e
%                    c
%                    s
%                     Data.(opt.epochArray{e}).(conditionArray{c}).(stopOutcomeArray{s}).sdf(i,:) = ...
%                         Data.(opt.epochArray{e}).(conditionArray{c}).(stopOutcomeArray{s}).sdf(i,:) / iNormFactor;
%                 end
%                 if ~strcmp(opt.epochArray{e}, 'responseOnset')
%
%                     Data.(opt.epochArray{e}).(conditionArray{c}).stopStop.sdf(i,:) = ...
%                         Data.(opt.epochArray{e}).(conditionArray{c}).stopStop.sdf(i,:) / iNormFactor;
%                 end
%             end
%         end
%     end
% end




%   ____________________ SET UP PLOT  ____________________

lineWidth = 3;  % for all conditions right now

nCol = length(opt.epochArray);
if ~opt.doStops
    nRow = 2;
else
    nRow = 2;
end

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
axStopStop = 1;
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
    
    if opt.doStops
        % Stop trials
        ax(axStop, e) = axes('units', 'centimeters', 'position', [xAxesPosition(axStop, e) yAxesPosition(axStop, e) axisWidth axisHeight]);
        hold(ax(axStop, e), 'on')
        %             set(ax(axStop, e), 'ylim', [yLimMin yLimMax], 'xlim', [epochRange(1) epochRange(end)])
        set(ax(axStop, e), 'ylim', [yLimMin yLimMax], 'xlim', [1 epochRange(end) - epochRange(1)])
        cla
        hold(ax(axStop, e), 'on')
        plot(ax(axStop, e), [-epochRange(1) -epochRange(1)], [yLimMin yLimMax * .9], '-k', 'linewidth', 2)
        set(ax(axStop, e), 'xticklabel', get(gca, 'xtick')+epochRange(1)-1)
        if e > 1, set(ax(axStop, e), 'yticklabel', []), end
        
        
        % Stop Stop trials
        ax(axStopStop, e) = axes('units', 'centimeters', 'position', [xAxesPosition(axStopStop, e) yAxesPosition(axStopStop, e) axisWidth axisHeight]);
        hold(ax(axStopStop, e), 'on')
        %             set(ax(axStopStop, e), 'ylim', [yLimMin yLimMax], 'xlim', [epochRange(1) epochRange(end)])
        set(ax(axStopStop, e), 'ylim', [yLimMin yLimMax], 'xlim', [1 epochRange(end) - epochRange(1)])
        cla
        hold(ax(axStopStop, e), 'on')
        plot(ax(axStopStop, e), [-epochRange(1) -epochRange(1)], [yLimMin yLimMax * .9], '-k', 'linewidth', 2)
        set(ax(axStopStop, e), 'xticklabel', get(gca, 'xtick')+epochRange(1)-1)
        if e > 1, set(ax(axStopStop, e), 'yticklabel', []), end
    end
end





% Loop through epochs
for e = 1 : length(opt.epochArray)
    epochRange = ccm_epoch_range(opt.epochArray{e}, 'plot');
    
    % Loop colorherence (easyIn, easyOut, hardIn, hardOut)
    for c = 1 : length(conditionArrayInd)
        cInd = conditionArrayInd(c);
        
        if opt.doStops
            %   _______ STOP TRIALS  _______
            % Index the relelvant units for the StopTarg Data
            stopTargInd = ismember(cell2table(Data.(opt.epochArray{e}).(conditionArray{cInd}).stopTarg.unitStop,'VariableNames',{'sessionID' 'unit'}), neurons(:,1:2));
            meanSDF = nanmean(Data.(opt.epochArray{e}).(conditionArray{cInd}).stopTarg.sdf(stopTargInd,:), 1);
            align = Data.(opt.epochArray{e}).alignStop;
            plot(ax(axStop, e), meanSDF(align+epochRange(1):end), 'LineWidth',lineWidth,'LineStyle',inOutStyle{cInd},'color',stopInOutColor(cInd,:))
            
            %         if ~strcmp(opt.epochArray{e}, 'stopSignalOn')
            %             meanSDF = nanmean(Data.(conditionArray{cInd}).goFast.(opt.epochArray{e}).sdf, 1);
            %             align = Data.(conditionArray{cInd}).goFast.(opt.epochArray{e}).align;
            %             plot(ax(axStop, e), meanSDF(align+epochRange(1):end), 'LineWidth',lineWidth,'LineStyle',inOutStyle{cInd},'color',goInOutColor(cInd,:))
            %         end
            
            
            %   _______ STOP STOP TRIALS  _______
            if ~strcmp(opt.epochArray{e}, 'responseOnset')
                % Index the relelvant units for the StopStop Data
                stopStopInd = ismember(cell2table(Data.(opt.epochArray{e}).(conditionArray{cInd}).stopStop.unitStop,'VariableNames',{'sessionID' 'unit'}), neurons(:,1:2));
                meanSDF = nanmean(Data.(opt.epochArray{e}).(conditionArray{cInd}).stopStop.sdf(stopStopInd,:), 1);
                align = Data.(opt.epochArray{e}).alignStop;
                plot(ax(axStopStop, e), meanSDF(align+epochRange(1):end), 'LineWidth',lineWidth,'LineStyle',inOutStyle{cInd},'color',stopStopColor(cInd,:))
            end
            %         if ~strcmp(opt.epochArray{e}, 'stopSignalOn')
            %             meanSDF = nanmean(Data.(conditionArray{cInd}).goSlow.(opt.epochArray{e}).sdf, 1);
            %             align = Data.(conditionArray{cInd}).goSlow.(opt.epochArray{e}).align;
            %             plot(ax(axStopStop, e), meanSDF(align+epochRange(1):end), 'LineWidth',lineWidth,'LineStyle',inOutStyle{cInd},'color',goInOutColor(cInd,:))
            %         end
        end  % opt.doStops

        
        
        %   _______ GO TRIALS  _______
        if ~strcmp(opt.epochArray{e}, 'stopSignalOn')
            % Go Targ
            meanSDF = nanmean(Data.(opt.epochArray{e}).(conditionArray{cInd}).goTarg.sdf(goTargInd,:), 1);
            %             align = Data.alignGo(e);
            align = Data.(opt.epochArray{e}).alignGo;
            plot(ax(axGo, e), meanSDF(align+epochRange(1):end), 'LineWidth',lineWidth,'LineStyle',inOutStyle{cInd},'color',goInOutColor(cInd,:))
            
            %                 meanSDF = nanmean(Data.(conditionArray{cInd}).goDist.(opt.epochArray{e}).sdf, 1);
            %                 align = Data.(conditionArray{cInd}).goDist.(opt.epochArray{e}).align;
            %                 plot(ax(axGo, e), meanSDF(align+epochRange(1):end), 'LineWidth',lineWidth,'LineStyle',inOutStyle{cInd},'color',goInOutColor(cInd,:))
        end
        
        
        
    end % colorCohArray
    
end % epochs


h=axes('Position', [0 0 1 1], 'Visible', 'Off');
titleString = sprintf('%s\t n = %d', opt.categoryName , size(neurons, 1));
text(0.5,1, titleString, 'HorizontalAlignment','Center', 'VerticalAlignment','Top')

if opt.printPlot
    filePath = fullfile(projectRoot,'results',projectDate,subject);
    if opt.easyOnly
        fileName = ['pop_',opt.categoryName,'_easy',addMulti,addNorm,'.pdf'];
    elseif ~opt.doStops
        fileName = ['pop_',opt.categoryName,'_go',addMulti,addNorm,'.pdf'];
    else
        fileName = ['pop_',opt.categoryName,addMulti,addNorm,'.pdf'];
    end
    print(figureHandle, fullfile(filePath, fileName), '-dpdf', '-r300')
end
return

% ________________________________________________________________________________
% ________________________________________________________________________________
% ________________________________________________________________________________
% ________________________________________________________________________________
% ________________________________________________________________________________



