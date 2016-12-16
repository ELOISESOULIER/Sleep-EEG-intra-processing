%compute average time frequency map over cycle

addpath(genpath('/home/eloise/Stages/Rythm-ICM/Matlab/Toolbox/'));
addpath(genpath('/home/eloise/Stages/Rythm-ICM/Matlab/cycle'));
addpath(genpath('/home/eloise/Stages/Rythm-ICM/Matlab/Code Matlab'));

params = getParams();
v_FreqAxis = linspace(100,params.s_freqMin,params.s_TFreso);
set1 = [{'1615'},{'1520'},{'2051'}];
set2 = [{'976'},{'1176'},{'1036'},{'1277'},{'2098'}];
mean1Evok = [];
mean1Spont = [];
mean2Evok = [];
mean2Spont = [];

for s_cell = 1:length(set1)
    
    load([set1{s_cell},'-Evok.mat'])
    load([set1{s_cell},'-Spont.mat'])
    
    if isempty(mean1Evok)
        mean1Evok = m_sumTFEvok;
        mean1Spont = m_sumTFSpont;
    else
        mean1Evok = mean1Evok + m_sumTFEvok;
        mean1Spont = mean1Spont + m_sumTFSpont;
    end
    
end

mean1Evok = mean1Evok./length(set1);
mean1Spont = mean1Spont./length(set1);

for s_cell = 1:length(set2)
    
    load([set2{s_cell},'-Evok.mat'])
    load([set2{s_cell},'-Spont.mat'])
    
    if isempty(mean2Evok)
        mean2Evok = m_sumTFEvok;
        mean2Spont = m_sumTFSpont;
    else
        mean2Evok = mean2Evok + m_sumTFEvok;
        mean2Spont = mean2Spont + m_sumTFSpont;
    end
    
end
mean2Evok = mean2Evok./length(set2);
mean2Spont = mean2Spont./length(set2);

v_TimeAxis = 1:0.36:length(mean1Evok)*0.36;

figure(1)
subplot(1,2,1)
    v_clim = [prctile(mean1Spont(:),1),prctile(mean1Spont(:),99)];
    f_ImageArray(mean1Spont,v_TimeAxis,v_FreqAxis,'colormap','jet','limits',v_clim)
    title('Spontané')
    colorbar
subplot(1,2,2)
    v_clim = [prctile(mean1Evok(:),1),prctile(mean1Evok(:),99)];
    f_ImageArray(mean1Evok,v_TimeAxis,v_FreqAxis,'colormap','jet','limits',v_clim)
    title('Évoqué')
    colorbar  

figure(2)
subplot(1,2,1)
    v_clim = [prctile(mean2Spont(:),1),prctile(mean2Spont(:),99)];
    f_ImageArray(mean2Spont,v_TimeAxis,v_FreqAxis,'colormap','jet','limits',v_clim)
    title('Spontané')
    colorbar
subplot(1,2,2)
    v_clim = [prctile(mean2Evok(:),1),prctile(mean2Evok(:),99)];
    f_ImageArray(mean2Evok,v_TimeAxis,v_FreqAxis,'colormap','jet','limits',v_clim)
    title('Évoqué')
    colorbar  