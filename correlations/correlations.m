%computes correlations between intra and EEG on spindle and 3 gamma
%frequency intervals : Low (30-60), Middle (60-80) and High(80-100)

info = openSpike2('period5.smr');
[m_upstates,number,info] = upstate_analysis(info);
s_Si = info(1).header.sampleinterval;
s_Fe = info(1).header.Sampling;

m_corrs = repmat(struct('number',0,'length',0,'corr_sp',0,'timelag_sp',0,...
    'corr_gamLo',0,'timelag_gamLo',0,'corr_gamMid',0,'timelag_gamMid',0,...
    'corr_gamHi',0,'timelag_gamHi',0),length(m_upstates),1);

for s_which=1:length(m_upstates)

%% Select "good" upstates and filter %%
    v_Vm = m_upstates{s_which,1};
    v_EEG = m_upstates{s_which,2};
    %test length and absence of action potentials
    b_valid = validity_test(v_Vm,s_Si);
    if b_valid
        %filter the 4 interesting frequency intervals
        [v_FiltSpVm,v_FiltGamLoVm,v_FiltGamMidVm, v_FiltGamHiVm,...
        v_FiltSpEEG,v_FiltGamLoEEG,v_FiltGamMidEEG,v_FiltGamHiEEG] =...
        filtre_upstate(v_Vm,v_EEG,info(1).header.Sampling);
        assert((length(v_FiltSpVm)==length(v_FiltGamLoVm)) && ...
                (length(v_FiltSpVm)==length(v_FiltGamMidVm)) && ...
                (length(v_FiltSpVm)==length(v_FiltGamHiVm)) && ...
                (length(v_FiltSpVm)==length(v_FiltSpEEG)) && ...
                (length(v_FiltSpEEG)==length(v_FiltGamLoEEG)) && ...
                (length(v_FiltSpEEG)==length(v_FiltGamMidEEG)) && ...
                (length(v_FiltSpEEG)==length(v_FiltGamHiEEG)));
        s_len = length(v_FiltSpVm);
        m_corrs(s_which).number = s_which;
        m_corrs(s_which).length = s_len*0.36;
        
%% visualize %%
        v_t = (0: length(v_FiltSpVm)-1)/s_Fe;
        figure(2)
        subplot(2,2,1);
            plot(v_t,v_FiltSpVm);
            title('Spindle 9 -20 vm');
            xlabel('Time (second-s)'); ylabel('mV');
        subplot(2,2,2);
            plot(v_t,v_FiltGamLoVm);
            title('Gamma 30 -60 vm');
            xlabel('Time (seconds)'); ylabel('mV');
        subplot(2,2,3);
            plot(v_t,v_FiltGamMidVm);
            title('Gamma 60 - 80 vm');
            xlabel('Time (seconds)'); ylabel('mV');
        subplot(2,2,4);
            plot(v_t,v_FiltGamHiVm);
            title('Gamma 80 - 100 vm');
            xlabel('Time (seconds)'); ylabel('mV');

        figure(3)
        subplot(2,2,1);
            plot(v_t,v_FiltSpEEG);
            title('Spindle 9 -20 EEG');
            xlabel('Time (seconds)'); ylabel('microV');
        subplot(2,2,2);
            plot(v_t,v_FiltGamLoEEG);
            title('Gamma 30 -60 EEG');
            xlabel('Time (seconds)'); ylabel('microV');
        subplot(2,2,3);
            plot(v_t,v_FiltGamMidEEG);
            title('Gamma 60 - 80 EEG');
            xlabel('Time (seconds)'); ylabel('microV');
        subplot(2,2,4);
            plot(v_t,v_FiltGamHiEEG);
            title('Gamma 80 - 100 EEG');
            xlabel('Time (seconds)'); ylabel('microV');

%% compute correlations and time lags %%

        x_sp = v_FiltSpVm - mean(v_FiltSpVm);
        y_sp = v_FiltSpEEG - mean(v_FiltSpEEG);
        x_gamLo = v_FiltGamLoVm - mean(v_FiltGamLoVm);
        y_gamLo = v_FiltGamLoEEG - mean(v_FiltGamLoEEG);
        x_gamMid = v_FiltGamMidVm - mean(v_FiltGamMidVm);
        y_gamMid = v_FiltGamMidEEG - mean(v_FiltGamMidEEG);
        x_gamHi = v_FiltGamHiVm - mean(v_FiltGamHiVm);
        y_gamHi = v_FiltGamHiEEG - mean(v_FiltGamHiEEG);

        [v_xcorrSp,v_lagSp] = xcorr(x_sp,y_sp, 'coeff');
        [v_xcorrGamLo,v_lagGamLo] = xcorr(x_gamLo,y_gamLo, 'coeff');
        [v_xcorrGamMid,v_lagGamMid] = xcorr(x_gamMid,y_gamMid, 'coeff');
        [v_xcorrGamHi,v_lagGamHi] = xcorr(x_gamHi,y_gamHi, 'coeff');
      
        [s_corrSp,I_sp] = max(abs(v_xcorrSp));
        s_timeDiffSp = v_lagSp(I_sp)/s_Fe;
        m_corrs(s_which).corr_sp = s_corrSp;
        m_corrs(s_which).timelag_sp = s_timeDiffSp;

        [s_corrGamLo,I_gamLo] = max(abs(v_xcorrGamLo));
        s_timeDiffGamLo = v_lagGamLo(I_gamLo)/s_Fe;
        m_corrs(s_which).corr_gamLo = s_corrGamLo;
        m_corrs(s_which).timelag_gamLo = s_timeDiffGamLo;

        [s_corrGamMid,I_gamMid] = max(abs(v_xcorrGamMid));
        s_timeDiffGamMid = v_lagGamMid(I_gamMid)/s_Fe;
        m_corrs(s_which).corr_gamMid = s_corrGamMid;
        m_corrs(s_which).timelag_gamMid = s_timeDiffGamMid;

        [s_corrGamHi,I_gamHi] = max(abs(v_xcorrGamHi));
        s_timeDiffGamHi = v_lagGamHi(I_gamHi)/s_Fe;
        m_corrs(s_which).corr_gamHi = s_corrGamHi;
        m_corrs(s_which).timelag_gamHi = s_timeDiffGamHi;
        
%% visualize %%
        figure(6)
        subplot(2,2,1);
        plot(v_lagSp,v_xcorrSp);
        title(strcat('Spindle : timediff =',num2str(s_timeDiffSp)));
        subplot(2,2,2);
        plot(v_lagGamLo, v_xcorrGamLo);
        title(strcat('Gamma 30-60 : timediff =',num2str(s_timeDiffGamLo)));
        subplot(2,2,3);
        plot(v_lagGamMid, v_xcorrGamMid);
        title(strcat('Gamma 60-80 : timediff =',num2str(s_timeDiffGamMid)));
        subplot(2,2,4);
        plot(v_lagGamHi, v_xcorrGamHi);
        title(strcat('Gamma 80-100 : timediff = ',num2str(s_timeDiffGamHi)));
    else continue
    end
end

