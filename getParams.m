%function intended to homogenize parameter values accross scripts

function params = getParams()

s_baseline = 2000; %segment before stim time in ms
s_TFreso = 200;  %Time frequency map resolution
s_nbBins = 400; %Number of bins in relative map
s_freqMin = 5;
s_freqMax = 100;
s_trialLen = 4000; %after stim time, in ms
params = struct('s_baseline',s_baseline,'s_TFreso',s_TFreso,...
    's_nbBins',s_nbBins,'s_freqMin',s_freqMin,'s_trialLen',s_trialLen,...
    's_freqMax',s_freqMax);

end
