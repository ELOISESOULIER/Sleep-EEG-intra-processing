%computes correlations between envelopes of intra and EEG signal on spindle
%and 3 gamma frequency intervals : Low (30-60), Middle (60-80) and High(80-100)

info = openSpike2('period7.smr');
[m_upstates,number,info] = upstate_analysis(info);
s_Si = info(1).header.sampleinterval;
s_Fe = info(1).header.Sampling;

m_corrsEnv = repmat(struct('number',0,'length',0,'corr_sp',0,'timelag_sp',0,...
    'corr_gamLo',0,'timelag_gamLo',0,'corr_gamMid',0,'timelag_gamMid',0,...
'corr_gamHi',0,'timelag_gamHi',0),length(m_upstates),1);

for s_which=1:length(m_upstates)
    %% Select "good" upstates and filter %%

    v_Vm= m_upstates{s_which,1};
    v_EEG = m_upstates{s_which,2};
    %test minimum length and absence of action potentials
    b_valid = validity_test(v_Vm,s_Si);
    if b_valid
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
        m_corrsEnv(s_which).number = s_which;
        m_corrsEnv(s_which).length = s_len*0.36;
        
%% compute signal envelope %%
        v_EnvSpVm = abs(hilbert(v_FiltSpVm'));
        v_EnvSpEEG = abs(hilbert(v_FiltSpEEG'));
        v_EnvGamLoVm = abs(hilbert(v_FiltGamLoVm'));
        v_EnvGamLoEEG = abs(hilbert(v_FiltGamLoEEG'));
        v_EnvGamMidVm = abs(hilbert(v_FiltGamMidVm'));
        v_EnvGamMidEEG = abs(hilbert(v_FiltGamMidEEG'));
        v_EnvGamHiVm = abs(hilbert(v_FiltGamHiVm'));
        v_EnvGamHiEEG = abs(hilbert(v_FiltGamHiEEG'));
        
%% visualize %%
        v_t = 0:0.36:0.36*(length(v_Vm));
        figure(11)
        subplot(2,2,1);
            plot(v_t,v_FiltSpVm);
            hold on;
            plot(v_t,[-1;1]*v_EnvSpVm,'r','LineWidth',2);
            title('Spindle 9 -20 vm');
            xlabel('Time (seconds)'); ylabel('mV');
            hold off;
        subplot(2,2,2);
            plot(v_t,v_FiltGamLoVm);
            hold on;
            plot(v_t,[-1;1]*v_EnvGamLoVm,'r','LineWidth',2);
            hold off;
            title('Gamma 30 -60 vm');
            xlabel('Time (seconds)'); ylabel('mV');
        subplot(2,2,3);
            plot(v_t,v_FiltGamMidVm);
            hold on;
            plot(v_t,[-1;1]*v_EnvGamMidVm,'r','LineWidth',2);
            hold off;
            title('Gamma 60 - 80 vm');
            xlabel('Time (seconds)'); ylabel('mV');
        subplot(2,2,4);
            plot(v_t,v_FiltGamHiVm);
            hold on;
            plot(v_t,[-1;1]*v_EnvGamHiVm,'r','LineWidth',2);
            hold off;
            title('Gamma 80 - 100 vm');
            xlabel('Time (seconds)'); ylabel('mV');

        figure(3)
        subplot(2,2,1);
            plot(v_t,v_FiltSpEEG);
            hold on;
            plot(v_t,[-1;1]*v_EnvSpEEG,'r','LineWidth',2);
            hold off;
            title('Spindle 9 -20 EEG');
            xlabel('Time (seconds)'); ylabel('microV');
        subplot(2,2,2);
            plot(v_t,v_FiltGamLoEEG);
            hold on;
            plot(v_t,[-1;1]*v_EnvGamLoEEG,'r','LineWidth',2);
            hold off;
            title('Gamma 30 -60 EEG');
            xlabel('Time (seconds)'); ylabel('microV');
        subplot(2,2,3);
            plot(v_t,v_FiltGamMidEEG);
            hold on;
            plot(v_t,[-1;1]*v_EnvGamMidEEG,'r','LineWidth',2);
            hold off;
            title('Gamma 60 - 80 EEG');
            xlabel('Time (seconds)'); ylabel('microV');
        subplot(2,2,4);
            plot(v_t,v_FiltGamHiEEG);
            hold on;
            plot(v_t,[-1;1]*v_EnvGamHiEEG,'r','LineWidth',2);
            hold off;
            title('Gamma 80 - 100 EEG');
            xlabel('Time (seconds)'); ylabel('microV');

%% compute correlations and time lag %%
        x_sp = v_EnvSpVm - mean(v_EnvSpVm);
        y_sp = v_EnvSpEEG - mean(v_EnvSpEEG);
        x_gamLo = v_EnvGamLoVm - mean(v_EnvGamLoVm);
        y_gamLo = v_EnvGamLoEEG - mean(v_EnvGamLoEEG);
        x_gamMid = v_EnvGamMidVm - mean(v_EnvGamMidVm);
        y_gamMid = v_EnvGamMidEEG - mean(v_EnvGamMidEEG);
        x_gamHi = v_EnvGamHiVm - mean(v_EnvGamHiVm);
        y_gamHi = v_EnvGamHiEEG - mean(v_EnvGamHiEEG);

        [v_xcorrSp,v_lagSp] = xcorr(x_sp,y_sp, 'coeff');
        [v_xcorrGamLo,v_lagGamLo] = xcorr(x_gamLo,y_gamLo, 'coeff');
        [v_xcorrGamMid,v_lagGamMid] = xcorr(x_gamMid,y_gamMid, 'coeff');
        [v_xcorrGamHi,v_lagGamHi] = xcorr(x_gamHi,y_gamHi, 'coeff');
      
        [s_corrSp,I_sp] = max(abs(v_xcorrSp));
        s_timeDiffSp = v_lagSp(I_sp)/s_Fe;
        m_corrsEnv(s_which).corr_sp = s_corrSp;
        m_corrsEnv(s_which).timelag_sp = s_timeDiffSp;

        [s_corrGamLo,I_gamLo] = max(abs(v_xcorrGamLo));
        s_timeDiffGamLo = v_lagGamLo(I_sp)/s_Fe;
        m_corrsEnv(s_which).corr_gamLo = s_corrGamLo;
        m_corrsEnv(s_which).timelag_gamLo = s_timeDiffGamLo;

        [s_corrGamMid,I_gamMid] = max(abs(v_xcorrGamMid));
        s_timeDiffGamMid = v_lagGamMid(I_gamMid)/s_Fe;
        m_corrsEnv(s_which).corr_gamMid = s_corrGamMid;
        m_corrsEnv(s_which).timelag_gamMid = s_timeDiffGamMid;

        [s_corrGamHi,I_gamHi] = max(abs(v_xcorrGamHi));
        s_timeDiffGamHi = v_lagGamHi(I_gamHi)/s_Fe;
        m_corrsEnv(s_which).corr_gamHi = s_corrGamHi;
        m_corrsEnv(s_which).timelag_gamHi = s_timeDiffGamHi;
        
%% Visualize %%
        
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
        %keyboard;
       
    end
end

 