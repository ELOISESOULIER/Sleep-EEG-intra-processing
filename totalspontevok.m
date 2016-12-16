% Computes all mean time frequency maps for one cell

addpath(genpath('/home/eloise/Stages/Rythm-ICM/Spike2-data/donnees_stimulees/'))
addpath(genpath('/home/eloise/Stages/Rythm-ICM/Matlab/Toolbox/'))

b_divided = true; %are spontaneous and evoked channels together 

if b_divided
    st_infoSpont = openSpike2('c1176.3-spontfqcy.smr');
    st_infoEvok = openSpike2('c1176.3-evokfqcy1.smr');
    s_sampling = st_infoSpont(2).header.Sampling;
    s_sampleInt = st_infoSpont(2).header.sampleinterval;
else
    st_info = openSpike2('c1307-fqcy.smr');
    s_sampling = st_info(2).header.Sampling;
    s_sampleInt = st_info(2).header.sampleinterval;
end

params = getParams();
s_TFreso = params.s_TFreso;
s_nbBins = params.s_nbBins;
s_freqMin = params.s_freqMin;
s_baseline = params.s_baseline;
v_FreqAxis = linspace(100,params.s_freqMin,s_TFreso);
v_FreqAxisTot = linspace(100,8,s_TFreso);


%% Time-Frequency maps
if b_divided %compute everything separately for evoked and spontaneous conditions
    [m_TFAbsVmSpont,m_TFAbsEEGSpont,~,~,v_TimeAxis]=...
         TFabs1sec_spontevok(st_infoSpont,true,false,s_TFreso,false);
    [m_TFRelVmSpont,m_TFRelEEGSpont,~,~]=...
        TFrelatif_spontevok(st_infoSpont,true,false,s_nbBins,s_TFreso,false);
    
    [~,~,m_TFAbsVmEvok,m_TFAbsEEGEvok,v_TimeAxis]=...
         TFabs1sec_spontevok(st_infoEvok,false,true,s_TFreso,false);
    [~,~,m_TFRelVmEvok,m_TFRelEEGEvok]=...
        TFrelatif_spontevok(st_infoEvok,false,true,s_nbBins,s_TFreso,false);
    [m_TFAbsVmSpontNorm,m_TFAbsEEGSpontNorm,~,~] = ...
        TFabs1sec_spontevok(st_infoSpont,true,false);
    [m_TFRelVmSpontNorm,m_TFRelEEGSpontNorm,~,~]=...
        TFrelatif_spontevok(st_infoSpont,true,false);  

    [~,~,m_TFAbsVmEvokNorm,m_TFAbsEEGEvokNorm] = ...
        TFabs1sec_spontevok(st_infoEvok,false,true);
    [~,~,m_TFRelVmEvokNorm,m_TFRelEEGEvokNorm]=...
        TFrelatif_spontevok(st_infoEvok,false,true);   
% 
%     [m_spont,~,st_infoSpont] = spont_evok_extract(st_infoSpont,true,false,true);
%     [~,m_evok,st_infoEvok] = spont_evok_extract(st_infoEvok,false,true,true);
%     s_sampInt = st_infoSpont(1).header.sampleinterval;

else
    [m_TFAbsVmSpont,m_TFAbsEEGSpont,m_TFAbsVmEvok,m_TFAbsEEGEvok,v_TimeAxis]=...
         TFabs1sec_spontevok(st_info,true,true,s_TFreso,false);
    [m_TFRelVmSpont,m_TFRelEEGSpont,m_TFRelVmEvok,m_TFRelEEGEvok]=...
        TFrelatif_spontevok(st_info,true,true,s_nbBins,s_TFreso,false);

    [m_TFAbsVmSpontNorm,m_TFAbsEEGSpontNorm,m_TFAbsVmEvokNorm,...
        m_TFAbsEEGEvokNorm]= TFabs1sec_spontevok(st_info);   
    [m_TFRelVmSpontNorm,m_TFRelEEGSpontNorm,m_TFRelVmEvokNorm,m_TFRelEEGEvokNorm]=...
         TFrelatif_spontevok(st_info);   
% 
%     [m_spont,m_evok,st_info] = spont_evok_extract(st_info,true,true,true);
%     s_sampInt = st_info(1).header.sampleinterval;
end

%% cut unwanted part of the baseline

s_sampleInt = 360;
s_cut = floor(400/s_sampleInt*1000);
v_TimeAxis = v_TimeAxis(1:end-s_cut+1);

m_TFAbsVmSpontNorm = m_TFAbsVmSpontNorm(:,s_cut:end);
m_TFAbsVmEvokNorm = m_TFAbsVmEvokNorm(:,s_cut:end);

m_TFAbsEEGSpontNorm = m_TFAbsEEGSpontNorm(:,s_cut:end);
m_TFAbsEEGEvokNorm = m_TFAbsEEGEvokNorm(:,s_cut:end);

m_TFAbsVmSpont = m_TFAbsVmSpont(:,s_cut:end);
m_TFAbsVmEvok = m_TFAbsVmEvok(:,s_cut:end);

m_TFAbsEEGSpont = m_TFAbsEEGSpont(:,s_cut:end);
m_TFAbsEEGEvok = m_TFAbsEEGEvok(:,s_cut:end);


%% Visualize
figure(20)
subplot(1,2,1)
    v_clim = [prctile(m_TFAbsVmSpontNorm(:),5),prctile(m_TFAbsVmSpontNorm(:),99)];
    f_ImageArray(m_TFAbsVmSpontNorm,v_TimeAxis,v_FreqAxisTot,'colormap','jet','limits',v_clim);
    title('TF absolu - Intra - Spontané');xlabel('Temps (s)');ylabel('Frequence (Hz)')
subplot(1,2,2)
    v_clim = [prctile(m_TFAbsVmSpontNorm(:),5),prctile(m_TFAbsVmSpontNorm(:),99)];
    f_ImageArray(m_TFAbsVmEvokNorm,v_TimeAxis,v_FreqAxisTot,'colormap','jet','limits',v_clim);
    title('TF absolu - Intra - Évoqué');xlabel('Temps (s)');ylabel('Frequence (Hz)')
    
% subplot(2,2,3)
%      v_clim = [prctile(m_TFRelVmSpontNorm(:),1),prctile(m_TFRelVmSpontNorm(:),99)];
%      f_ImageArray(m_TFRelVmSpontNorm,(1:1:s_nbBins)*1/s_nbBins,v_FreqAxisTot,...
%          'colormap','jet','limits',v_clim);
%      title('TF relatif - Intra - Spontané');xlabel('Temps (%upstate)');ylabel('Frequence (Hz)')
% subplot(2,2,4)
%      v_clim = [prctile(m_TFRelVmEvokNorm(:),1),prctile(m_TFRelVmEvokNorm(:),99)];
%      f_ImageArray(m_TFRelVmEvokNorm,(1:1:s_nbBins)*1/s_nbBins,v_FreqAxisTot,...
%          'colormap','jet','limits',v_clim);
%      title('TF relatif - Intra - Évoqué');xlabel('Temps (%upstate)');ylabel('Frequence (Hz)')
    
figure(21)
subplot(1,2,1)
    v_clim = [prctile(m_TFAbsEEGSpontNorm(:),5),prctile(m_TFAbsEEGSpontNorm(:),99)];
    f_ImageArray(m_TFAbsEEGSpontNorm,v_TimeAxis,v_FreqAxisTot,'colormap','jet','limits',v_clim);
    title('TF absolu - EEG - Spontané');xlabel('Temps (s)');ylabel('Frequence (Hz)')
    save('/home/eloise/Stages/Rythm-ICM/Matlab/1176-trials/totalSpontNorm','m_TFAbsEEGSpontNorm');

subplot(1,2,2)
    v_clim = [prctile(m_TFAbsEEGEvokNorm(:),5),prctile(m_TFAbsEEGEvokNorm(:),99)];
    f_ImageArray(m_TFAbsEEGEvokNorm,v_TimeAxis,v_FreqAxisTot,'colormap','jet','limits',v_clim);
    title('TF absolu - EEG - Évoqué');xlabel('Temps (s)');ylabel('Frequence (Hz)')
    saveas(gcf,'/home/eloise/Stages/Rythm-ICM/Matlab/1176-trials/totalNorm.png');
    save('/home/eloise/Stages/Rythm-ICM/Matlab/1176-trials/totalEvokNorm','m_TFAbsEEGEvokNorm');

% subplot(2,2,3)
%     v_clim = [prctile(m_TFRelEEGSpontNorm(:),1),prctile(m_TFRelEEGSpontNorm(:),99)];
%     f_ImageArray(m_TFRelEEGSpontNorm,(1:1:s_nbBins)*1/s_nbBins,v_FreqAxisTot,...
%         'colormap','jet','limits',v_clim);
%     title('TF relatif - EEG - Spontané');xlabel('Temps (%upstate)');ylabel('Frequence (Hz)')
% subplot(2,2,4)
%     v_clim = [prctile(m_TFRelEEGEvokNorm(:),1),prctile(m_TFRelEEGEvokNorm(:),99)];
%     f_ImageArray(m_TFRelEEGEvokNorm,(1:1:s_nbBins)*1/s_nbBins,v_FreqAxisTot,...
%         'colormap','jet','limits',v_clim);
%     title('TF relatif - EEG - Évoqué');xlabel('Temps (%upstate)');ylabel('Frequence (Hz)')

figure(31)
subplot(1,2,1)
    v_clim = [prctile(m_TFAbsVmSpont(:),5),prctile(m_TFAbsVmSpont(:),99)];
    f_ImageArray(m_TFAbsVmSpont,v_TimeAxis,v_FreqAxis,'colormap','jet','limits',v_clim);
    title('TF absolu - Intra - Spontané');xlabel('Temps (s)');ylabel('Frequence (Hz)')
subplot(1,2,2)
    v_clim = [prctile(m_TFAbsVmEvok(:),5),prctile(m_TFAbsVmEvok(:),99)];
    f_ImageArray(m_TFAbsVmEvok,v_TimeAxis,v_FreqAxis,'colormap','jet','limits',v_clim);
    title('TF absolu - Intra - Évoqué');xlabel('Temps (s)');ylabel('Frequence (Hz)')
% subplot(2,2,3)
%      v_clim = [prctile(m_TFRelVmSpont(:),1),prctile(m_TFRelVmSpont(:),99)];
%      f_ImageArray(m_TFRelVmSpont,(1:1:s_nbBins)*1/s_nbBins,v_FreqAxis,...
%          'colormap','jet','limits',v_clim);
%      title('TF relatif - Intra - Spontané');xlabel('Temps (%upstate)');ylabel('Frequence (Hz)')
% subplot(2,2,4)
%      v_clim = [prctile(m_TFRelVmEvok(:),1),prctile(m_TFRelVmEvok(:),99)];
%      f_ImageArray(m_TFRelVmEvok,(1:1:s_nbBins)*1/s_nbBins,v_FreqAxis,...
%          'colormap','jet','limits',v_clim);
%      title('TF relatif - Intra - Évoqué');xlabel('Temps (%upstate)');ylabel('Frequence (Hz)') 
%     
figure(32)
subplot(1,2,1)
    v_clim = [prctile(m_TFAbsEEGSpont(:),5),prctile(m_TFAbsEEGSpont(:),99)];
    f_ImageArray(m_TFAbsEEGSpont,v_TimeAxis,v_FreqAxis,'colormap','jet','limits',v_clim);
    title('TF absolu - EEG - Spontané');xlabel('Temps (s)');ylabel('Frequence (Hz)')
    save('/home/eloise/Stages/Rythm-ICM/Matlab/1176-trials/totalSpont','m_TFAbsEEGSpont')

subplot(1,2,2)
    v_clim = [prctile(m_TFAbsEEGEvok(:),5),prctile(m_TFAbsEEGEvok(:),99)];
    f_ImageArray(m_TFAbsEEGEvok,v_TimeAxis,v_FreqAxis,'colormap','jet','limits',v_clim);
    title('TF absolu - EEG - Évoqué');xlabel('Temps (s)');ylabel('Frequence (Hz)')                     
    saveas(gcf,'/home/eloise/Stages/Rythm-ICM/Matlab/1176-trials/total.png');
    save('/home/eloise/Stages/Rythm-ICM/Matlab/1176-trials/totalEvok','m_TFAbsEEGEvok')

   
%% normalization factor for EEG
    s_meanEEG = 0;
    s_stdEEG = 0; 
    s_nbEEG = 0;
    for s_which = 1:length(m_spont)
        b_valid = true;%validity_test(m_spont{s_which,1},st_info(1).header.sampleinterval);
        if b_valid
            s_meanEEG = s_meanEEG + mean(m_spont{s_which,2});
            s_stdEEG = s_stdEEG + std(m_spont{s_which,2});
            s_nbEEG = s_nbEEG + 1;
        end
    end
    s_meanEEG = s_meanEEG/s_nbEEG;
    s_stdEEG = s_stdEEG/s_nbEEG;
    
%% Raw and filtered signal    
%     s_countSpont = 0;
% for s_which=1:length(m_spont)
%     %[b_valid,m_spont{s_which,1}] = validity_test(m_spont{s_which,1},...
%      %   st_info(1).header.sampleinterval);
%     b_valid = true;
%     if b_valid
%         s_countSpont = s_countSpont + 1;
%         m_spont{s_which,1} = m_spont{s_which,1} -70 - m_spont{s_which,1}(1);
%         m_spont{s_which,2} = (m_spont{s_which,2}-s_meanEEG)/s_stdEEG^2;
%         m_spont{s_which,1} = m_spont{s_which,1}(floor(400/s_sampleInt*1000):end);
%         m_spont{s_which,2} = m_spont{s_which,2}(floor(400/s_sampleInt*1000):end);
%         [v_FiltSpVm,v_FiltGamLoVm,v_FiltGamMidVm, v_FiltGamHiVm,...
%             v_FiltSpEEG,v_FiltGamLoEEG,v_FiltGamMidEEG,v_FiltGamHiEEG] =...
%             filtre_upstate(m_spont{s_which,1},m_spont{s_which,2},s_sampling);
%         figure(27)
%         time = 0:s_sampleInt/1000000:s_sampleInt/1000000*(length(m_spont{s_which,1})-1);
%         Ay1(1)=subplot(5,2,1);
%             plot(time, m_spont{s_which,1}); hold on;
%             title('Intra - Spontané'); ylabel('Potentiel (mV)')
%         Ay2(1)=subplot(5,2,3);
%             plot(time, v_FiltSpVm(1:end-1));hold on;
%             title('Intra - Spindle');ylabel('Potentiel(mV)')
%         Ay3(1)=subplot(5,2,5)
%             plot(time,v_FiltGamLoVm(1:end-1));hold on;
%             title('Intra - Gamma 30-60');ylabel('Potentiel(mV)')
%         Ay4(1)=subplot(5,2,7)
%             plot(time,v_FiltGamMidVm(1:end-1));hold on;
%             title('Intra - Gamma 60-80');ylabel('Potentiel(mV)')
%         Ay5(1)=subplot(5,2,9)
%             plot(time,v_FiltGamHiVm(1:end-1));hold on;
%             title('Intra - Gamma 80-100');ylabel('Potentiel(mV)')
%        figure(28)     
%        Ay6(1)=subplot(5,2,1);
%             plot(time, m_spont{s_which,2});hold on;
%             title('EEG - Spontané');ylabel('Potentiel (uV)')
%        Ay7(1)=subplot(5,2,3);
%             plot(time, v_FiltSpEEG(1:end-1));hold on;
%             title('EEG - Spindle');ylabel('Potentiel(uV)')
%        Ay8(1)=subplot(5,2,5);
%             plot(time,v_FiltGamLoEEG(1:end-1));hold on;
%             title('EEG - Gamma 30-60');ylabel('Potentiel (uV)')     
%        Ay9(1)=subplot(5,2,7);
%             plot(time,v_FiltGamMidEEG(1:end-1));hold on;
%             title('EEG - Gamma 60-80');ylabel('Potentiel(uV)')
%        Ay10(1)=subplot(5,2,9);
%             plot(time,v_FiltGamHiEEG(1:end-1));hold on;
%             title('EEG - Gamma 80-100');ylabel('Potentiel (uV)')
%             
%     end
% end
% 
% 
% s_countEvok = 0;    
% 
% for s_which=1:length(m_evok)
%     %[b_valid,m_evok{s_which,1}] = validity_test(m_evok{s_which,1},...
%      %   st_info(1).header.sampleinterval,true);
%     b_valid = true;
%     if b_valid
%         s_countEvok = s_countEvok + 1;
%         m_evok{s_which,1} = m_evok{s_which,1}-70 -m_evok{s_which,1}(1);
%         m_evok{s_which,2} = (m_evok{s_which,2}-s_meanEEG)/s_stdEEG^2;
%         m_evok{s_which,1} = m_evok{s_which,1}(floor(400/s_sampleInt*1000):end);
%         m_evok{s_which,2} = m_evok{s_which,2}(floor(400/s_sampleInt*1000):end);
%         [v_FiltSpVm,v_FiltGamLoVm,v_FiltGamMidVm, v_FiltGamHiVm,...
%             v_FiltSpEEG,v_FiltGamLoEEG,v_FiltGamMidEEG,v_FiltGamHiEEG] =...
%             filtre_upstate(m_evok{s_which,1},m_evok{s_which,2},s_sampling);
%         
%         time = 0:s_sampleInt/1000000:s_sampleInt/1000000*(length(m_evok{s_which,1})-1);
%         
%         figure(27)
%         Ay1(2)=subplot(5,2,2);
%             plot(time, m_evok{s_which,1});hold on;
%             title('Intra - Évoqué');ylabel('Potentiel (mV)')
%         Ay2(2)=subplot(5,2,4);
%             plot(time, v_FiltSpVm(1:end-1));hold on;
%             title('Intra - Spindle');ylabel('Potentiel(mV)')
%         Ay3(2)=subplot(5,2,6);
%             plot(time,v_FiltGamLoVm(1:end-1));hold on;
%             title('Intra - Gamma 30-60');ylabel('Potentiel(mV)')
%         Ay4(2)=subplot(5,2,8);
%             plot(time,v_FiltGamMidVm(1:end-1));hold on;
%             title('Intra - Gamma 60-80');ylabel('Potentiel(mV)')
%         Ay5(2)=subplot(5,2,10);
%             plot(time,v_FiltGamHiVm(1:end-1));hold on;
%             title('Intra - Gamma 80-100');ylabel('Potentiel(mV)')
%       figure(28)
%         Ay6(2)=subplot(5,2,2);
%             plot(time, m_evok{s_which,2});hold on;
%             title('EEG - Évoqué');ylabel('Potentiel (uV)')
%         Ay7(2)=subplot(5,2,4);
%             plot(time, v_FiltSpEEG(1:end-1));hold on;
%             title('EEG - Spindle');ylabel('Potentiel(uV)')
%         Ay8(2)=subplot(5,2,6);
%             plot(time,v_FiltGamLoEEG(1:end-1));hold on;
%             title('EEG - Gamma 30-60');ylabel('Potentiel (uV)')
%         Ay9(2)=subplot(5,2,8);
%             plot(time,v_FiltGamMidEEG(1:end-1));hold on;
%             title('EEG - Gamma 60-80');ylabel('Potentiel(uV)')
%         Ay10(2)=subplot(5,2,10);
%             plot(time,v_FiltGamHiEEG(1:end-1));hold on;
%             title('EEG - Gamma 80-100');ylabel('Potentiel (uV)')
%         linkaxes(Ay1,'y'); linkaxes(Ay2,'y'); linkaxes(Ay3,'y');
%         linkaxes(Ay4,'y'); linkaxes(Ay5,'y'); linkaxes(Ay6,'y');
%         linkaxes(Ay7,'y'); linkaxes(Ay8,'y'); linkaxes(Ay9,'y');
%         linkaxes(Ay10,'y');
%         
%     end
% end
