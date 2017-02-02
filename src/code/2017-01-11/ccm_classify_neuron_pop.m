function neuronTypes = ccm_classify_neuron_pop(subject,projectRoot,projectDate, append)
%
% Create a table (using stats from all sessions) of sessions with neurons, classifying the neurons w.r.t different epochs:
%
%
dataPath = fullfile(projectRoot,'data',projectDate,subject);


% Load the list of units that was created using ccm_list_units.m
load(fullfile(dataPath, 'ccm_units.mat'))

sessionID = units.sessionID;
unit = units.unit;
 
if append
    load(fullfile(dataPath, 'ccm_neuronTypes'), 'neuronTypes')
    neuronTypes = table2cell(neuronTypes);
    
%     lastSession = neuronTypes(end, 1);
    %     ind = find(strcmp(lastSession, neuronTypes.sessionID));
%     startInd = 1 + find(strcmp(lastSession, sessionList));
    startInd = 1 + size(neuronTypes, 1);
    
    % Add empty cells the size of the new units to add
    neuronTypes = [neuronTypes; cell(length(sessionID) - size(neuronTypes, 1), 13)];
else
 	neuronTypes = cell(size(units, 1), 13);
    startInd = 1;
end

opt                 = ccm_options;
opt.collapseSignal  = true;
opt.collapseTarget  = true;
opt.doStops         = false;
opt.plotFlag        = false;
opt.printPlot       = false;

tic
poolID = parpool(6);
parfor i = startInd :  length(sessionID)
% startInd = 298;
% for i = startInd :  startInd%length(sessionID)
    
    
    fprintf('%02d\t%s\t%s\n',i,sessionID{i}, unit{i})
    
    iUnit               = [sessionID(i), unit(i)];
    iData               = ccm_session_data(subject, iUnit, opt);
    iCat              	= ccm_classify_neuron(iData);
    neuronTypes(i,:)  	= iCat;
    

end
delete(poolID)
toc
neuronTypes = cell2table(neuronTypes, 'VariableNames',...
    {'sessionID'...
    'unit'...
    'hemisphere'...
    'rf'...
    'fix'...
    'vis'...
    'checker'...
    'presacc'...
    'presaccMax'...
    'presaccRamp'...
    'postsacc'...
    'reward'...
    'intertrial'});

save(fullfile(dataPath, 'ccm_neuronTypes'), 'neuronTypes')


