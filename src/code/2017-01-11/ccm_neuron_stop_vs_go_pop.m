function ccm_neuron_stop_vs_go_pop(subject,projectRoot,projectDate, options)
%
% Create a table (using stats from all sessions) of sessions with neurons, classifying the neurons w.r.t different epochs:
%
% dataPath = fullfile(projectRoot,'data',projectDate,subject);
%
%

ANALYZE_CANCELED = options.ANALYZE_CANCELED;
ANALYZE_NONCANCELED = options.ANALYZE_NONCANCELED;



if options.multiUnit
    addMulti = '_multiUnit';
else
    addMulti = [];
end
dataPath = fullfile(projectRoot,'data',projectDate,subject);
fileName = ['ccm_neuronTypes',addMulti];


original = load(fullfile(dataPath, fileName));

sessionID = original.neuronTypes.sessionID;
unit = original.neuronTypes.unit;

% testLength = 30;
% sessionID = sessionID(1:testLength);
% unit = unit(1:testLength);



if options.append
    %     % Load one of the population data files and determine last session
    %     % entered
    %     load(fullfile(dataPath, 'ccm_canceled_vs_go_neuronTypes'), 'cancelTypes')
    %     load(fullfile(dataPath, 'ccm_noncanceled_vs_go_neuronTypes'), 'noncancelTypes')
    %     lastSession = cancelTypes.sessionID(end);
    %     lastSession = noncancelTypes.sessionID(end);
    %     startInd = 1 + find(strcmp(lastSession, sessionList), 1);
    %
    % %         startInd = 1 + size(cancelGoNeuronData, 1);
    %
    % %     disp(cancelGoNeuronData(end-10:end,:))
    %     fprintf('\nAppending starting on next session\n')
else
    cancelTypes = cell(length(sessionID), 13);
    noncancelTypes =  cell(length(sessionID), 9);
    %     cancelTypes = cell(testLength, 12);
    %     noncancelTypes =  cell(testLength, 9);
    startInd = 1;
end
tic
poolID = parpool(options.parpoolSize);
parfor i = startInd : length(sessionID)
    % for i = startInd : length(sessionID)
    %     fprintf('%02d\t%s\t%s\n',i,sessionID{i}, unit{i})
    
    iData = ccm_neuron_stop_vs_go(subject, sessionID{i}, unit(i), options);
    
    
    if ANALYZE_CANCELED
        % Collect the Canceled vs Go trials
        iCCell = {sessionID(i), ...
            unit(i), ...
            {iData.rf}, ...
            iData.stopStopSsrt, ...
            iData.stopStopCoh, ...
            iData.stopStopSsd, ...
            iData.stopStopCond, ...
            iData.nStopStop, ...
            iData.pValue40msStopStop, ...
            iData.cancelTimeDist, ...
            iData.cancelTime2Std, ...
            iData.cancelTime4Std, ...
            iData.cancelTime6Std};
        
        cancelTypes(i,:) = iCCell;
    end
    
    
    
    if ANALYZE_NONCANCELED
        
        % Collect the Noncanceled vs Go trials
        iNcCell = {sessionID(i), ...
            unit(i), ...
            {iData.rf}, ...
            iData.stopTargSsrt, ...
            iData.stopTargCoh, ...
            iData.stopTargSsd, ...
            iData.stopTargCond, ...
            iData.nStopTarg, ...
            iData.pValue40msStopTarg};
        
        
        noncancelTypes(i,:) = iNcCell;
    end
    
    
end
delete(poolID)
toc


if ANALYZE_CANCELED
    cancelTypes = cell2table(cancelTypes, 'VariableNames',...
        {'sessionID'...
        'unit'...
        'rf'...
        'stopStopSsrt'...
        'stopStopCoh'...
        'stopStopSsd'...
        'stopStopCond'...
        'nStopStop'...
        'pValue40msStopStop'...
        'cancelTimeDist'...
        'cancelTime2Std'...
        'cancelTime4Std'...
        'cancelTime6Std'});
    
    if ~isdir(fullfile(dataPath, 'go_vs_canceled', options.ssrt))
        mkdir(fullfile(dataPath, 'go_vs_canceled', options.ssrt))
    end
    save(fullfile(dataPath, 'go_vs_canceled', options.ssrt, ['ccm_canceled_vs_go_neuronTypes',addMulti]), 'cancelTypes')
end


if ANALYZE_NONCANCELED
    noncancelTypes = cell2table(noncancelTypes, 'VariableNames',...
        {'sessionID'...
        'unit'...
        'rf'...
        'stopTargSsrt'...
        'stopTargCoh'...
        'stopTargSsd'...
        'stopTargCond'...
        'nStopTarg'...
        'pValue40msStopTarg'});
    
    if ~isdir(fullfile(dataPath, 'go_vs_noncanceled', options.ssrt))
        mkdir(fullfile(dataPath, 'go_vs_noncanceled', options.ssrt))
    end
    save(fullfile(dataPath, 'go_vs_noncanceled', options.ssrt, ['ccm_noncanceled_vs_go_neuronTypes',addMulti]), 'noncancelTypes')
end

