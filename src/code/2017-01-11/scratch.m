
%%
cancelData.cancelTime = cancelData.cancelTime2Std - cancelData.stopStopSsd - cancelData.ssrt;

test = cellfun(@(x,y,z) x - y - z, cancelData.cancelTime2Std, cancelData.stopStopSsd, num2cell(cancelData.ssrt), 'uni', false)
