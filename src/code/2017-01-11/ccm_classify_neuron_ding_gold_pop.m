function ccm_classify_neuron_ding_gold_pop(subject,projectRoot,projectDate, append)
%
% Create a table (using stats from all sessions) of sessions with neurons, classifying the neurons w.r.t different epochs:
% Stim, Sacc, Post, Reward
%
%
dataPath = fullfile(projectRoot,'data',projectDate,subject);
original = load(fullfile(dataPath, 'ccm_neuronTypes'));

sessionID = original.neuronTypes.sessionID;
unit = original.neuronTypes.unit;
rf = original.neuronTypes.rf;

% Either append the data to an extant file, or create a new table and start
% from the beginning
if append
    %     load(fullfile(dataPath, 'ccm_ding_gold_neuronTypes'), 'neuronTypes')
    %
    %     % Figure out last unit being processed/saved
    %     lastSession = neuronTypes.sessionID(end);
    %     lastUnit = neuronTypes.unit(end);
    %     startInd = 1 + find(strcmp(original.neuronTypes.sesssionID, lastSession) & strcmp(original.neuronTypes.unit, lastUnit));
else
    stimTypes = cell(size(unit, 1), 12);
    saccTypes = cell(size(unit, 1), 12);
    postTypes = cell(size(unit, 1), 12);
    rewardTypes = cell(size(unit, 1), 12);
    startInd = 1;
    %     startUnit = 1;
end


opt                 = ccm_options;
opt.doStops         = false;
opt.plotFlag        = false;
opt.printPlot       = false;


% Loop through the sessions and add the data to the table.
poolID = parpool(6);
parfor i = startInd : size(original.neuronTypes, 1)
% for i = startInd : size(original.neuronTypes, 1)
    
    
        fprintf('%02d\t%s\t%s\n',i,sessionID{i}, unit{i})
    
    iUnit = [sessionID(i), unit(i)];
    iData               = ccm_session_data(subject, iUnit, opt);
    iData.rf       = rf{i};
    
    iData.epoch = 'Stim';
    iStim           = ccm_classify_neuron_ding_gold(iData);
    stimTypes(i,:) = iStim;
    
    iData.epoch = 'Sacc';
    iSacc           = ccm_classify_neuron_ding_gold(iData);
    saccTypes(i,:) = iSacc;
    
    iData.epoch = 'Post';
    iPost           = ccm_classify_neuron_ding_gold(iData);
    postTypes(i,:) = iPost;
    
%     iData.epoch = 'Reward';
%     iRew           = ccm_classify_neuron_ding_gold(iData);
%     rewardTypes(i,:) = iRew;
    
end
delete(poolID)

neuronTypes = cell2table(stimTypes, 'VariableNames',...
    {'sessionID'...
    'unit'...
    'hemisphere'...
    'rf'...
    'epoch'...
    'choice'...
    'coherence'...
    'ddm'...
    'tChoice'...
    'leftIsIn'...
    'coeffIn'...
    'coeffOut'});

save(fullfile(dataPath, 'ccm_ddmStim_neuronTypes'), 'neuronTypes')

neuronTypes = cell2table(saccTypes, 'VariableNames',...
    {'sessionID'...
    'unit'...
    'hemisphere'...
    'rf'...
    'epoch'...
    'choice'...
    'coherence'...
    'ddm'...
    'tChoice'...
    'leftIsIn'...
    'coeffIn'...
    'coeffOut'});

save(fullfile(dataPath, 'ccm_ddmSacc_neuronTypes'), 'neuronTypes')

neuronTypes = cell2table(postTypes, 'VariableNames',...
    {'sessionID'...
    'unit'...
    'hemisphere'...
    'rf'...
    'epoch'...
    'choice'...
    'coherence'...
    'ddm'...
    'tChoice'...
    'leftIsIn'...
    'coeffIn'...
    'coeffOut'});

save(fullfile(dataPath, 'ccm_ddmPost_neuronTypes'), 'neuronTypes')

% neuronTypes = cell2table(rewardTypes, 'VariableNames',...
%     {'sessionID'...
%     'unit'...
%     'hemisphere'...
%     'rf'...
%     'epoch'...
%     'choice'...
%     'coherence'...
%     'ddm'...
%     'tChoice'...
%     'leftIsIn'...
%     'coeffIn'...
%     'coeffOut'});
% 
% save(fullfile(dataPath, 'ccm_ddmReward_neuronTypes'), 'neuronTypes')
