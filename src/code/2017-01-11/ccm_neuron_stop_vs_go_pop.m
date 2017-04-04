function ccm_neuron_stop_vs_go_pop(subject,projectRoot,projectDate, options)
%
% Create a table (using stats from all sessions) of sessions with neurons, classifying the neurons w.r.t different epochs:
%
% dataPath = fullfile(projectRoot,'data',projectDate,subject);
%
%

ANALYZE_NONCANCELED = false;



dataPath = fullfile(projectRoot,'data',projectDate,subject);
fileName = 'ccm_neuronTypes';
if options.multiUnit
    fileName = [fileName, '_multiUnit'];
end

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
    cancelTypes = cell(length(sessionID), 12);
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
        iData.cancelTime2Std, ...
        iData.cancelTime4Std, ...
        iData.cancelTime6Std};
    
    cancelTypes(i,:) = iCCell;
    
    
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
    'cancelTime2Std'...
    'cancelTime4Std'...
    'cancelTime6Std'});

if ~isdir(fullfile(dataPath, 'go_vs_canceled', options.ssrt))
    mkdir(fullfile(dataPath, 'go_vs_canceled', options.ssrt))
end
if options.multiUnit
    save(fullfile(dataPath, 'go_vs_canceled', options.ssrt, 'ccm_canceled_vs_go_neuronTypes_multiUnit'), 'cancelTypes')
else
    save(fullfile(dataPath, 'go_vs_canceled', options.ssrt, 'ccm_canceled_vs_go_neuronTypes'), 'cancelTypes')
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
    if options.multiUnit
        save(fullfile(dataPath, 'go_vs_noncanceled', options.ssrt, 'ccm_noncanceled_vs_go_neuronTypes_multiUnit'), 'noncancelTypes')
    else
        save(fullfile(dataPath, 'go_vs_noncanceled', options.ssrt, 'ccm_noncanceled_vs_go_neuronTypes'), 'noncancelTypes')
    end
end

