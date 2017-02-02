function ccm_classify_neuron_ding_gold_pop(subject,projectRoot,projectDate, append)
%
% Create a table (using stats from all sessions) of sessions with neurons, classifying the neurons w.r.t different epochs:
% Stim, Sacc, Post, Reward
%
%
dataPath = fullfile(projectRoot,'data',projectDate,subject);


% Open the sessions file and makes lists of the entries
fid=  fopen(fullfile(dataPath,['ccm_sessions_',subject,'.csv']));


nCol = 5;
formatSpec = '%s';
mHeader = textscan(fid, formatSpec, nCol, 'Delimiter', ',');

mData = textscan(fid, '%s %s %d %d %d', 'Delimiter', ',','TreatAsEmpty',{'NA','na'});

sessionList     = mData{1};
hemisphereList  = mData{2};
neuronLogical   = mData{3};



% Determine options for calling ccm_session_data.m
opt             = ccm_options;
opt.howProcess  = 'each';
opt.plotFlag    = false;
opt.printPlot    = false;
opt.dataType    = 'neuron';
opt.collapseTarg 	= true;
opt.collapseSignal 	= false;
opt.doStops 	= false;




% Either append the data to an extant file, or create a new table and start
% from the beginning
if append
    load(fullfile(dataPath, 'ccm_ding_gold_neuronTypes'), 'neuronTypes')
    
    % Figure out last session being processed/saved
    lastSession = neuronTypes.sessionID(end);
    startInd = 1 + find(strcmp(lastSession, sessionList));

    % Figure out last unit saved 
    lastUnit = neuronTypes.unit(end);
    sessionInd =find(strcmp(lastSession, neuronTypes.sessionID));
    unitSess = unique(neuronTypes.unit(sessionInd));
    startUnit = 1 + find(strcmp(neuronTypes.unit(end), unitSess));
else
    neuronTypes = table();
    startInd = 1;
    startUnit = 1;
end





% Loop through the sessions and add the data to the table.
for i = startInd : length(sessionList)
    
    
    if neuronLogical(i)
        fprintf('\n%02d\t%s\n',i,sessionList{i})
        
        % See how many units we'll loop through for this session (to save
        % disk space  - so matlab doesn't crash)
        [td, S, E] = load_data(subject, sessionList{i});
        nUnit = length(S.spikeUnitArray);
        
        opt.hemisphere = hemisphereList{i};
        opt.trialData = td;
        opt.SessionData = S;
        opt.ExtraVar = E;
        
       
        for j = startUnit : nUnit
            fprintf('\t%02d\t%s\n',j,S.spikeUnitArray{j})
            opt.unitArray   = S.spikeUnitArray(j);
            iData           = ccm_session_data(subject, sessionList{i}, opt);
            iData.hemisphere = opt.hemisphere;
            jUnit           = ccm_classify_neuron_ding_gold(iData);
            
            neuronTypes     = [neuronTypes; jUnit];
            clear iData
       end
            save(fullfile(dataPath, 'ccm_ding_gold_neuronTypes'), 'neuronTypes')
    end
end
