function ccm_neuron_stop_vs_go_pop(subject,projectRoot,projectDate, append)
%
% Create a table (using stats from all sessions) of sessions with neurons, classifying the neurons w.r.t different epochs:
%
dataPath = fullfile(projectRoot,'data',projectDate,subject);


% Open the table of neurons classified
allTypes = load(fullfile(dataPath, 'ccm_neuronTypes'));
allTypes = allTypes.neuronTypes;


opt             = ccm_options;
opt.howProcess  = 'each';
opt.plotFlag    = false;
opt.printPlot    = false;
opt.dataType    = 'neuron';
opt.collapseTarg 	= true;
opt.minTrialPerCond 	= 10;

sessionList = unique(allTypes.sessionID);

if append
    % Load one of the population data files and determine last session
    % entered
%     load(fullfile(dataPath, 'ccm_canceled_vs_go_neuronTypes'), 'cancelGoNeuronData')
    load(fullfile(dataPath, 'ccm_noncanceled_vs_go_neuronTypes'), 'noncancelGoNeuronData')
%     lastSession = cancelGoNeuronData.sessionID(end);
    lastSession = noncancelGoNeuronData.sessionID(end);
    startInd = 1 + find(strcmp(lastSession, sessionList), 1);
    
%         startInd = 1 + size(cancelGoNeuronData, 1);

%     disp(cancelGoNeuronData(end-10:end,:))
%     fprintf('\nAppending starting on next session\n')
else
    cancelGoNeuronData = table();
%     noncancelGoNeuronData = table();
    startInd = 1;
end

for i = startInd : length(sessionList)
    fprintf('%02d\t%s\n',i,sessionList{i})
  
    % See how many units we'll loop through for this session (to save
    % disk space  - so matlab doesn't crash)
    [~, S, ~] = ccm_load_data_behavior(subject, sessionList{i});
    nUnit = length(S.spikeUnitArray);
    
    
    
    for j = 1 : nUnit
       fprintf('\t%02d\t%s\n',j,S.spikeUnitArray{j})
        
        iData = ccm_neuron_stop_vs_go(subject, sessionList{i}, S.spikeUnitArray(j), opt);
        
        % Collect the Canceled vs Go trials
        cUnit = table();
        cUnit.sessionID     = repmat(sessionList(i), length(iData.stopStopSsd), 1);
        cUnit.rf          = repmat(iData.rf, length(iData.stopStopSsd), 1);
        cUnit.ssrt          = repmat(iData.ssrt, length(iData.stopStopSsd), 1);
        cUnit.unit          = repmat(S.spikeUnitArray(j), length(iData.stopStopSsd), 1);
        cUnit.stopStopCoh    = iData.stopStopCoh;
        cUnit.stopStopSsd    = iData.stopStopSsd;
        
        cUnit.pValue40msStopStop    = iData.pValue40msStopStop;
        cUnit.cancelTime2Std    = iData.cancelTime2Std;
        cUnit.cancelTime4Std    = iData.cancelTime4Std;
        cUnit.cancelTime6Std    = iData.cancelTime6Std;
               
        cancelGoNeuronData = [cancelGoNeuronData; cUnit];

        
%         % Collect the Noncanceled vs Go trials
%         ncUnit = table();
%         ncUnit.sessionID     = repmat(sessionList(i), length(iData.stopTargSsd), 1);
%         ncUnit.hemisphere          = repmat(iData.rf, length(iData.hemisphere), 1);
%         ncUnit.rf          = repmat(iData.rf, length(iData.stopTargSsd), 1);
%         ncUnit.ssrt          = repmat(iData.ssrt, length(iData.stopTargSsd), 1);
%         ncUnit.unit          = repmat(S.spikeUnitArray(j), length(iData.stopTargSsd), 1);
%         ncUnit.stopTargCoh    = iData.stopTargCoh;
%         ncUnit.stopTargSsd    = iData.stopTargSsd;
%         
%         ncUnit.pValue40ms    = iData.pValue40msStopTarg;
%         
%         
%         noncancelGoNeuronData = [noncancelGoNeuronData; ncUnit];

neuronTypes = cancelGoNeuronData;
save(fullfile(dataPath, 'ccm_canceled_vs_go_neuronTypes'), 'neuronTypes')
% neuronTypes = noncancelGoNeuronData;
%         save(fullfile(dataPath, 'ccm_noncanceled_vs_go_neuronTypes'), 'neuronTypes')

        
        clear iData
    end
    
end
