function Data = ccm_population_neuron(subjectID,projectRoot,projectDate,Opt)

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
task = 'ccm';

dataPath = fullfile(projectRoot,'data',projectDate,subjectID);

if nargin < 2
    Opt.dataType        = 'neuron';
    Opt.sessionSet     	= 'neural1';
    
    Opt.sessionArray   	= [];
    Opt.unitArray    	= [];
    Opt.rfList          = [];
    Opt.hemisphereList 	= [];
    
    Opt.doStops       	= true;
    Opt.multiUnit       	= false;
    Opt.normalize       	= false;
    
    if nargin == 0
        Data = Opt;
        return
    end
end
MIN_TRIAL = 10;

sessionSet      = Opt.sessionSet;
sessionArray    = Opt.sessionArray;
dataType        = Opt.dataType;
rfList        = Opt.rfList;
hemisphereList        = Opt.hemisphereList;

unitArray = Opt.unitArray;


if isempty(sessionArray)
    [sessionArray, ~] = task_session_array(subjectID, task, sessionSet);
end

if Opt.multiUnit
    addMulti = '_multiUnit';
else
    addMulti = [];
end
if Opt.normalize
    addNorm = '_normalized';
else
    addNorm = [];
end
iNormFactor = 1;


nUnit = length(sessionArray);





% Define epochs, outcome categories, and color coherence conditions to
% collect
epochArrayStop      = {'fixWindowEntered', 'targOn', 'checkerOn', 'stopSignalOn', 'responseOnset', 'toneOn', 'rewardOn'};
epochArrayGo        = {'fixWindowEntered', 'targOn', 'checkerOn', 'responseOnset', 'toneOn', 'rewardOn'};
outcomeArrayGo      = {'goTarg', 'goDist'};
outcomeArrayStop    = {'stopTarg', 'stopStop','goFast', 'goSlow'};
colorCohArray       = {'easyIn', 'easyOut', 'hardIn', 'hardOut'};

nSsd = 20; % Make sure we have enough ssd slots per unit

unit    = cell(nUnit, 2);

rtGo    = nan(nUnit, length(colorCohArray), length(outcomeArrayGo));
rtStop  = nan(nUnit, length(colorCohArray), length(outcomeArrayGo), nSsd);

sdfGo    = cell(nUnit, length(epochArrayGo), length(colorCohArray), length(outcomeArrayGo));
sdfStop  = cell(nUnit, length(epochArrayStop), length(colorCohArray), length(outcomeArrayStop), nSsd);

ssdStop = nan(nUnit, length(epochArrayStop), length(colorCohArray), length(outcomeArrayStop), nSsd);



% Figure out which SSDs to collapse for go/stop comparison:
% Implement latency matching in ccm_session_data and send it here

%     % Intialize Data structure for go trials
for k = 1 : length(colorCohArray)
    for m = 1 : length(epochArrayGo)
        Data.(epochArrayGo{m}).(colorCohArray{k}).goTarg.sdf = [];
        Data.(epochArrayGo{m}).(colorCohArray{k}).goDist.sdf = [];
    end
    Data.checkerOn.(colorCohArray{k}).goTarg.rt = [];
    Data.checkerOn.(colorCohArray{k}).goDist.rt = [];
    %     % Intialize Data structure for stop trials
    for m = 1 : length(epochArrayStop)
        Data.(epochArrayStop{m}).(colorCohArray{k}).stopTarg.unitStop = [];
        Data.(epochArrayStop{m}).(colorCohArray{k}).stopDist.unitStop = [];
        Data.(epochArrayStop{m}).(colorCohArray{k}).stopStop.unitStop = [];
        Data.(epochArrayStop{m}).(colorCohArray{k}).goFast.unitStop = [];
        Data.(epochArrayStop{m}).(colorCohArray{k}).goSlow.unitStop = [];
        
        Data.(epochArrayStop{m}).(colorCohArray{k}).stopTarg.sdf = [];
        Data.(epochArrayStop{m}).(colorCohArray{k}).stopDist.sdf = [];
        Data.(epochArrayStop{m}).(colorCohArray{k}).stopStop.sdf = [];
        Data.(epochArrayStop{m}).(colorCohArray{k}).goFast.sdf = [];
        Data.(epochArrayStop{m}).(colorCohArray{k}).goSlow.sdf = [];
        
        Data.(epochArrayStop{m}).(colorCohArray{k}).stopTarg.ssd = [];
        Data.(epochArrayStop{m}).(colorCohArray{k}).stopDist.ssd = [];
        Data.(epochArrayStop{m}).(colorCohArray{k}).stopStop.ssd = [];
        Data.(epochArrayStop{m}).(colorCohArray{k}).goFast.ssd = [];
        Data.(epochArrayStop{m}).(colorCohArray{k}).goSlow.ssd = [];
    end
    Data.checkerOn.(colorCohArray{k}).stopTarg.rt = [];
    Data.checkerOn.(colorCohArray{k}).stopDist.rt = [];
    Data.checkerOn.(colorCohArray{k}).goFast.rt = [];
    Data.checkerOn.(colorCohArray{k}).goSlow.rt = [];
end







for m = 1 : length(epochArrayGo)
    mEpoch = epochArrayGo{m};
    sdfWindow = ccm_epoch_range(mEpoch, 'analyze');
    Data.(epochArrayGo{m}).alignGo = -sdfWindow(1);
    
end
for m = 1 : length(epochArrayStop)
    mEpoch = epochArrayStop{m};
    sdfWindow = ccm_epoch_range(mEpoch, 'analyze');
    Data.(epochArrayStop{m}).alignStop = -sdfWindow(1);
end





for iUnit = 1 : nUnit
    
    % Print out progress of the function
    fprintf('%d of %d\tUnit: %s\t%s\n', iUnit, nUnit, sessionArray{iUnit}, unitArray{iUnit})
    
    % Get single session neural data for this unit
    optData = ccm_options;
    optData.multiUnit = Opt.multiUnit;
    optData.printPlot = false;
    optData.plotFlag = false;
    sessionUnit = [sessionArray(iUnit), unitArray(iUnit)];
    iData = ccm_session_data(subjectID, sessionUnit, optData);
    
    ssdArray = iData.ssdArray;
    
    
    % Get the Receptive Field. If there's no RF, use the contralateral
    % direction relative to the recorded hemisphere.
    if strcmp(rfList{iUnit}, 'none')
        switch lower(hemisphereList{iUnit})
            case 'left'
                rfList{iUnit} = 'right';
            case 'right'
                rfList{iUnit} = 'left';
        end
    end
    
    
    % Figure out the indices of left and right color coherence proportions
    % for choices of hardest and easiest, into and out of RF (easyIn,
    % easyOut, hardIn, hardOut).
    pSignalArray = iData(1).pSignalArray;
    pSignalArray(pSignalArray == .5) = [];
    
    switch lower(rfList{iUnit})
        case 'left'
            easyInInd     = 1;
            easyOutInd    = length(pSignalArray);
            hardInInd     = length(pSignalArray)/2;
            hardOutInd    = hardInInd + 1;
        case 'right'
            easyInInd     = length(pSignalArray);
            easyOutInd    = 1;
            hardOutInd    = length(pSignalArray)/2;
            hardInInd     = hardOutInd + 1;
        case 'none'
            
    end
    iRFList = [easyInInd, easyOutInd, hardInInd, hardOutInd]; % Make sure this same order as colorCohArray
    
             if Opt.normalize
%         iNormFactor = max(iData.responseOnset.colorCoh(iRFList(1)).goTarg.sdfMean);
        normalWindow = -200 : 0;
        normalAlign = iData.targOn.colorCoh(iRFList(1)).goTarg.alignTime;
        iNormFactor = mean(iData.targOn.colorCoh(iRFList(1)).goTarg.sdfMean(normalAlign+normalWindow));
            end
    
    
    
    
    
    switch dataType
        case 'neuron'
            
            
            for k = 1 : length(colorCohArray)
                
               
                
                
                %              Collect Go Data
                % ============================================
                for n = 1 : length(outcomeArrayGo)
                    nOutcome = outcomeArrayGo{n};
                    
                    
                    % Collect mean RTs for the outcome
                    rtGo(iUnit, k, n) = ...
                        nanmean(iData.checkerOn.colorCoh(iRFList(k)).(nOutcome).rt);
                    
                    for m = 1 : length(epochArrayGo)
                        mEpoch = epochArrayGo{m};
                        sdfWindow = ccm_epoch_range(mEpoch, 'analyze');
                        
                        % Fill data with NaN if the condition has no trials
                        if size(iData.(mEpoch).colorCoh(iRFList(k)).(nOutcome).sdf, 1) <= MIN_TRIAL
                            sdfGo{iUnit, m, k, n} = ...
                                nan(1, 1+sdfWindow(end)-sdfWindow(1));
                        else
                            alignTime = iData.(mEpoch).colorCoh(iRFList(k)).(nOutcome).alignTime;
                            
                            % Might need to pad the sdf if aligntime is before the sdf window beginning
                            if alignTime <  -sdfWindow(1)
                                iData.(mEpoch).colorCoh(iRFList(k)).(nOutcome).sdf = ...
                                    [nan(size(iData.(mEpoch).colorCoh(iRFList(k)).(nOutcome).sdf, 1), -sdfWindow(1) - alignTime),...
                                    iData.(mEpoch).colorCoh(iRFList(k)).(nOutcome).sdf];
                                alignTime = alignTime - sdfWindow(1);
                            end
                            % Might need to pad the sdf if aligntime
                            % doesn't leave enough space before end of sdf
                            if alignTime + sdfWindow(end) > size(iData.(mEpoch).colorCoh(iRFList(k)).(nOutcome).sdf, 2)
                                iData.(mEpoch).colorCoh(iRFList(k)).(nOutcome).sdf = ...
                                    [iData.(mEpoch).colorCoh(iRFList(k)).(nOutcome).sdf,...
                                    nan(size(iData.(mEpoch).colorCoh(iRFList(k)).(nOutcome).sdf, 1), sdfWindow(end) - sdfWindow(1) - size(iData.(mEpoch).colorCoh(iRFList(k)).(nOutcome).sdf, 2)+1)];
                            end
                            sdfGo{iUnit, m, k, n} = ...
                                nanmean(iData.(mEpoch).colorCoh(iRFList(k)).(nOutcome).sdf(:,alignTime + sdfWindow), 1) ./ iNormFactor;
                        end
                    end
                    
                end
                
                
                
                
                
                
                
                %              Collect Stop Data
                % ============================================
                if Opt.doStops
                    for n = 1 : length(outcomeArrayStop)
                        nOutcome = outcomeArrayStop{n};
                        
                        for m = 1 : length(epochArrayStop)
                            mEpoch = epochArrayStop{m};
                            
                            % Skip collecting stop signal aligned data for
                            % latency-matched go trials
                            if ~((strcmp(nOutcome, 'goFast') || strcmp(nOutcome, 'goSlow')) && strcmp(mEpoch, 'stopSignalOn'))
                                sdfWindow = ccm_epoch_range(mEpoch, 'analyze');
                                
                                for s = 1 : length(ssdArray)
                                    % Are there enough trials in this stop
                                    % outcome?
                                    
                                    switch nOutcome
                                        case {'stopTarg', 'goFast'}
                                            sTrial = size(iData.(mEpoch).colorCoh(iRFList(k)).stopTarg.ssd(s).sdf, 1);
                                            gTrial = size(iData.(mEpoch).colorCoh(iRFList(k)).goFast.ssd(s).sdf, 1);
                                        case {'stopStop', 'goSlow'}
                                            if ~strcmp(mEpoch, 'responseOnset') % No RT data on stopStop trials
                                                sTrial = size(iData.(mEpoch).colorCoh(iRFList(k)).stopStop.ssd(s).sdf, 1);
                                                gTrial = size(iData.(mEpoch).colorCoh(iRFList(k)).goSlow.ssd(s).sdf, 1);
                                            else
                                            end
                                            %                                         case {'stopDist'}
                                            %                                             sTrial = size(iData.(mEpoch).colorCoh(iRFList(k)).stopDist.ssd(s).sdf, 1);
                                    end
                                    if sTrial >= MIN_TRIAL
                                        % Skip data collection if trying to collect
                                        % saccade-aligned data for stopStop trials
                                        if ~(strcmp(mEpoch, 'responseOnset') && strcmp(nOutcome, 'stopStop'));
                                            
                                            % Collect RTs if it's checker Onset epoch
                                            % and it's not stopStop outcome (make sure
                                            % saccade was made)
                                            if strcmp(mEpoch, 'checkerOn') && ~(strcmp(nOutcome, 'stopStop') || strcmp(nOutcome, 'goSlow'))
                                                rtStop(iUnit, k, n, s) = ...
                                                    nanmean(iData.(mEpoch).colorCoh(iRFList(k)).(nOutcome).ssd(s).rt);
                                            end
                                            
                                            alignTime = iData.(mEpoch).colorCoh(iRFList(k)).(nOutcome).ssd(s).alignTime;
                                            if ~isempty(alignTime)
                                                %                                 alignTime = nan;
                                                %                             end
                                                
                                                % Might need to pad the sdf if aligntime is before the sdf window beginning
                                                if alignTime <  -sdfWindow(1)
                                                    iData.(mEpoch).colorCoh(iRFList(k)).stopTarg.ssd(s).sdf = ...
                                                        [nan(size(iData.(mEpoch).colorCoh(iRFList(k)).(nOutcome).ssd(s).sdf, 1), -sdfWindow(1)+1 - alignTime),...
                                                        iData.(mEpoch).colorCoh(iRFList(k)).(nOutcome).ssd(s).sdf];
                                                    alignTime = alignTime - sdfWindow(1)+1;
                                                end
                                                % Might need to pad the sdf if aligntime
                                                % doesn't leave enough space before end of sdf
                                                if alignTime + sdfWindow(end) > size(iData.(mEpoch).colorCoh(iRFList(k)).(nOutcome).ssd(s).sdf, 2)
                                                    iData.(mEpoch).colorCoh(iRFList(k)).(nOutcome).ssd(s).sdf = ...
                                                        [iData.(mEpoch).colorCoh(iRFList(k)).(nOutcome).ssd(s).sdf,...
                                                        nan(size(iData.(mEpoch).colorCoh(iRFList(k)).(nOutcome).ssd(s).sdf, 1), sdfWindow(end) -sdfWindow(1) - size(iData.(mEpoch).colorCoh(iRFList(k)).(nOutcome).ssd(s).sdf, 2)+1)];
                                                end
                                                
                                                ssdStop(iUnit, m, k, n, s) = ...
                                                    ssdArray(s);
                                                %                                             mEpoch
                                                %                                             nOutcome
                                                %                                             s
                                                %                                             alignTime
                                                % if isempty(alignTime)
                                                %     disp('paused here')
                                                % end
                                                sdfStop{iUnit, m, k, n, s} = ...
                                                    nanmean(iData.(mEpoch).colorCoh(iRFList(k)).(nOutcome).ssd(s).sdf(:,alignTime + sdfWindow), 1) ./ iNormFactor;
                                            end
                                        end
                                    end
                                end
                            end
                            
                            
                        end
                    end
                end % if opt.doStops
                
                
            end % for k = 1 : length(colorCohArray)
            
        case 'lfp'
            
        case 'erp'
    end % swtich dataType
    
    
end % for iUnit = 1 : nSession



Data.unitGo = [sessionArray, unitArray];


% GO Data
for k = 1 : length(colorCohArray)
    for n = 1 : length(outcomeArrayGo)
        for m = 1 : length(epochArrayGo)
            
            Data.(epochArrayGo{m}).(colorCohArray{k}).(outcomeArrayGo{n}).sdf = cell2mat(sdfGo(:, m, k, n));
        end
        Data.checkerOn.(colorCohArray{k}).(outcomeArrayGo{n}).rt = rtGo(:, k, n);
    end
end

% STOP Data
Data.unitStop = [];
for iUnit = 1 : size(unit, 1)
    for k = 1 : length(colorCohArray)
        for n = 1 : length(outcomeArrayStop)
            nOutcome = outcomeArrayStop{n};
            for s = 1 : nSsd
                if ~isempty(sdfStop{iUnit, 1, k, n, s})
                    %                     Data.unitStop = [Data.unitStop; [sessionArray(iUnit), unitArray(iUnit)]];
                    for m = 1 : length(epochArrayStop)
                        mEpoch = epochArrayStop{m};
                        Data.(epochArrayStop{m}).(colorCohArray{k}).(outcomeArrayStop{n}).unitStop = ...
                            [Data.(epochArrayStop{m}).(colorCohArray{k}).(outcomeArrayStop{n}).unitStop; ...
                            [sessionArray(iUnit), unitArray(iUnit)]];
                        Data.(epochArrayStop{m}).(colorCohArray{k}).(outcomeArrayStop{n}).sdf = ...
                            [Data.(epochArrayStop{m}).(colorCohArray{k}).(outcomeArrayStop{n}).sdf; ...
                            sdfStop{iUnit, m, k, n, s}];
                        Data.(epochArrayStop{m}).(colorCohArray{k}).(outcomeArrayStop{n}).ssd = ...
                            [Data.(epochArrayStop{m}).(colorCohArray{k}).(outcomeArrayStop{n}).ssd; ...
                            ssdStop(iUnit, m, k, n, s)];
                        if strcmp(mEpoch, 'checkerOn') && ~(strcmp(nOutcome, 'stopStop') || strcmp(nOutcome, 'goSlow'))
                            Data.checkerOn.(colorCohArray{k}).(outcomeArrayStop{n}).rt = ...
                                [Data.checkerOn.(colorCohArray{k}).(outcomeArrayStop{n}).rt;...
                                rtStop(:, k, n, s)];
                        end
                    end
                end
            end
        end
    end
end




save(fullfile(dataPath, ['ccm_',Opt.categoryName,'_neuron_population',addMulti,addNorm]), '-struct', 'Data','-v7.3')
return


