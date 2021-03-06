function sessionList = ccm_exclude_sessions(subject)

switch subject
    case 'broca'
        sessionList = {...
            'bp080n01', ...
            'bp082n02', ...
            'bp084n02', ...
            'bp091n02-pm', ...
            'bp162n02'};
    case 'joule'
        sessionList = {...
            'jp054n02', ...
            'jp060n02', ...
            'jp061n02', ...
            'jp098n02', ...
            'jp104n02'};
end
            
            