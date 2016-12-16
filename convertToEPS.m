%open all .mat files containing time frequency maps and export them as eps
%images

addpath(genpath('C:\Users\Eloise\Documents\donnees_stimulees\'));
addpath(genpath('C:\Users\Eloise\Documents\Toolbox\'));

v_cells = [{'1176'};{'976'};{'1036'};{'1045'};{'1054'};{'1277'};...
    {'1307'};{'1520'};{'1615'};{'2051'};{'2098'};];

params = getParams();
v_TimeAxis = 0:0.36e-3:1.6;
v_FreqAxis = linspace(100,params.s_freqMin,params.s_TFreso);

for s_cellNum = 1:length(v_cells)    

    load(['C:\Users\Eloise\Documents\MATLAB\' v_cells{s_cellNum} '-trials/totalSpontNorm']);
    load(['C:\Users\Eloise\Documents\MATLAB\' v_cells{s_cellNum} '-trials/totalEvokNorm']);
    
    load(['C:\Users\Eloise\Documents\MATLAB\' v_cells{s_cellNum} '-trials/totalSpont']);
    load(['C:\Users\Eloise\Documents\MATLAB\' v_cells{s_cellNum} '-trials/totalEvok']);
    
    
     v_clim = [prctile(m_TFAbsEEGSpontNorm(:),5),prctile(m_TFAbsEEGSpontNorm(:),99)];
    figure(12)
    f_ImageArray(m_TFAbsEEGEvokNorm,v_TimeAxis,v_FreqAxis,'colormap','jet','limits',v_clim);
    title('TF normalis� - EEG - Evoqu�');xlabel('Temps (s)');ylabel('Frequence (Hz)')
    colorbar

    figure(13)
    f_ImageArray(m_TFAbsEEGSpontNorm,v_TimeAxis,v_FreqAxis,'colormap','jet','limits',v_clim);
    title('TF normalis� - EEG - Spontan�');xlabel('Temps (s)');ylabel('Frequence (Hz)')
    colorbar


    v_clim = [prctile(m_TFAbsEEGSpont(:),5),prctile(m_TFAbsEEGSpont(:),99)];
    figure(14)  
    f_ImageArray(m_TFAbsEEGEvok,v_TimeAxis,v_FreqAxis,'colormap','jet','limits',v_clim);
    title(' EEG - Evoqu�');xlabel('Temps (s)');ylabel('Frequence (Hz)')
    colorbar

    figure(15)
    f_ImageArray(m_TFAbsEEGSpont,v_TimeAxis,v_FreqAxis,'colormap','jet','limits',v_clim);
    title(' EEG - Spontan�');xlabel('Temps (s)');ylabel('Frequence (Hz)')
    colorbar

    %keyboard;
end
    