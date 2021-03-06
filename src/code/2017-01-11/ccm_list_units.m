function ccm_list_units(subject,projectRoot,projectDate, append, multiUnit)
%
% Create a table (using stats from all sessions) of sessions with neurons, classifying the neurons w.r.t different epochs:
%
% Need to create a csv file from manually entered sessions called
% ccm_sessions_<subject>.csv and save it into the appropriate data path for
% this date and subject.
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
lfpLogical      = mData{4};
eegLogical      = mData{5};


if append
    load(fullfile(dataPath, 'ccm_units'), 'units')
    lastSession = units.sessionID(end);
    
    startInd = 1 + find(strcmp(lastSession, sessionList));
else
    units = table();    
    %     neuronTypes = cell(length(sessionList), 13);
    startInd = 1;
end


for i = startInd :  length(sessionList)
    
    
    if neuronLogical(i)
        fprintf('%02d\t%s\n',i,sessionList{i})
        
        % See how many units we'll loop through for this session (to save
        % disk space  - so matlab doesn't crash)
        [~, S, ~] = load_data(subject, sessionList{i}, ccm_min_vars, multiUnit);
        jUnit = table();
        jUnit.sessionID = sessionList(i);
        jUnit.hemisphere = hemisphereList(i);
        for j = 1 : length(S.spikeUnitArray)
            jUnit.unit = S.spikeUnitArray(j);
            
            units = [units; jUnit];
        end
        if multiUnit
        save(fullfile(dataPath, 'ccm_multiunits'), 'units')
        else
        save(fullfile(dataPath, 'ccm_units'), 'units')
        end
    end
end
