function sessionKeep = ccm_session_by_electrode(subject,projectRoot,projectDate, electrode, channels)

dataPath = fullfile(projectRoot,'data',projectDate,subject);

number = table();

% Total channels recorded, total modulation neurons
load(fullfile(dataPath, 'ccm_neuronTypes'))
sessions = unique(neuronTypes.sessionID);
for i = 1 : length(sessions)
    iNumber = table();
    
    iNumber.session = sessions(i);
    % find the index of the first spike unit from this session
    firstInd = find(strcmp(sessions{i}, neuronTypes.sessionID), 1, 'first');
    lastInd = find(strcmp(sessions{i}, neuronTypes.sessionID), 1, 'last');
    
    % How many spike units recorded on each channel?
    unitInd = cellfun(@(x) str2double(regexp(x,'\d*','Match')), neuronTypes.unit(firstInd:lastInd));
    
    % How many channels/electrodes were recorded/used during this session?
    iChannels = unique(unitInd);
    iNumber.channels = length(iChannels);
    
    iNumber.unitPerChannel = cell(1,1);
    iNumber.unitPerChannel{1} = nan(1, iNumber.channels);
    
    for j = 1 : iNumber.channels
        iNumber.unitPerChannel{1}(j) = sum(iChannels(j) == unitInd);
    end
    
    iNumber.unitPerSession = sum(iNumber.unitPerChannel{1});
    number = [number; iNumber];
    
end


sessionKeep = number.session(ismember(number.channels, channels));

