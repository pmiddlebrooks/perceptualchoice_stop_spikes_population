function neuronTypes = ccm_classify_neuron_pop(subject,projectRoot,projectDate, options)
%
% Create a table (using stats from all sessions) of sessions with neurons, classifying the neurons w.r.t different epochs:
%
%
if nargin < 4
    options.multiUnit   = false;
    options.append      = false;
    options.parpoolSize = 4;
end

dataPath = fullfile(projectRoot,'data',projectDate,subject);


% Load the list of units that was created using ccm_list_units.m
switch options.multiUnit
    case false
        load(fullfile(dataPath, 'ccm_units.mat'))
    case true
        load(fullfile(dataPath, 'ccm_multiunits.mat'))
end

sessionID = units.sessionID;
unit = units.unit;

if options.append
    load(fullfile(dataPath, 'ccm_neuronTypes'), 'neuronTypes')
    neuronTypes = table2cell(neuronTypes);
    
    %     lastSession = neuronTypes(end, 1);
    %     ind = find(strcmp(lastSession, neuronTypes.sessionID));
    %     startInd = 1 + find(strcmp(lastSession, sessionList));
    startInd = 1 + size(neuronTypes, 1);
    
    % Add empty cells the size of the new units to add
    neuronTypes = [neuronTypes; cell(length(sessionID) - size(neuronTypes, 1), 14)];
else
    neuronTypes = cell(size(units, 1), 15);
    startInd = 1;
end






opt                 = ccm_options;
opt.collapseSignal  = true;
opt.collapseTarget  = true;
opt.doStops         = false;
opt.plotFlag        = true;
opt.printPlot       = true;
opt.figureHandle   	= 106;
opt.multiUnit       = options.multiUnit;




tic
poolID = parpool(options.parpoolSize);
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
    'saccadeBaseRatio'...
    'presaccMax'...
    'presaccRamp'...
    'presaccPeak'...
    'postsacc'...
    'reward'...
    'intertrial'});

if options.multiUnit
    save(fullfile(dataPath, 'ccm_neuronTypes_multiUnit'), 'neuronTypes')
else
    save(fullfile(dataPath, 'ccm_neuronTypes'), 'neuronTypes')
end



