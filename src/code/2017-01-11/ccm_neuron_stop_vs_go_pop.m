function ccm_neuron_stop_vs_go_pop(subject,projectRoot,projectDate, append)
%
% Create a table (using stats from all sessions) of sessions with neurons, classifying the neurons w.r.t different epochs:
%
% dataPath = fullfile(projectRoot,'data',projectDate,subject);
%
%





dataPath = fullfile(projectRoot,'data',projectDate,subject);
original = load(fullfile(dataPath, 'ccm_neuronTypes'));

sessionID = original.neuronTypes.sessionID;
unit = original.neuronTypes.unit;
rf = original.neuronTypes.rf;


opt             = ccm_options;
opt.howProcess  = 'each';
opt.plotFlag    = false;
opt.printPlot    = false;
opt.minTrialPerCond 	= 10;


if append
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
    cancelTypes = cell(length(sessionID), 10);
    noncancelTypes =  cell(length(sessionID), 7);
    startInd = 1;
end
tic
poolID = parpool(6);
parfor i = startInd : length(sessionID)
%     fprintf('%02d\t%s\t%s\n',i,sessionID{i}, unit{i})
    
    iData = ccm_neuron_stop_vs_go(subject, sessionID{i}, unit(i), opt);
    
    % Collect the Canceled vs Go trials
    iCCell = {sessionID(i), ...
        unit(i), ...
        {iData.rf}, ...
        iData.ssrt, ...
        iData.stopStopCoh, ...
        iData.stopStopSsd, ...
        iData.pValue40msStopStop, ...
        iData.cancelTime2Std, ...
        iData.cancelTime4Std, ...
        iData.cancelTime6Std};
    
    cancelTypes(i,:) = iCCell;
    
    
    
    % Collect the Noncanceled vs Go trials
    iNcCell = {sessionID(i), ...
        unit(i), ...
        {iData.rf}, ...
        iData.ssrt, ...
        iData.stopTargCoh, ...
        iData.stopTargSsd, ...
        iData.pValue40msStopTarg};
    
    
    noncancelTypes(i,:) = iNcCell;
    
    
    
end
delete(poolID)
toc

cancelTypes = cell2table(cancelTypes, 'VariableNames',...
    {'sessionID'...
    'unit'...
    'rf'...
    'ssrt'...
    'stopStopCoh'...
    'stopStopSsd'...
    'pValue40msStopStop'...
    'cancelTime2Std'...
    'cancelTime4Std'...
    'cancelTime6Std'});

noncancelTypes = cell2table(noncancelTypes, 'VariableNames',...
    {'sessionID'...
    'unit'...
    'rf'...
    'ssrt'...
    'stopTargCoh'...
    'stopTargSsd'...
    'pValue40msStopTarg'});

save(fullfile(dataPath, 'ccm_canceled_vs_go_neuronTypes'), 'cancelTypes')
save(fullfile(dataPath, 'ccm_noncanceled_vs_go_neuronTypes'), 'noncancelTypes')
