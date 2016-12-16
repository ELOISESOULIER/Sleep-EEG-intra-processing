%computes correlations between intra and EEG on 4 different frequency
%intervals for evoked and spontaneous signal

v_filesEvok = [{'c976.4-fqcy.smr'};{'c1036.2_fqcy.smr'};{'c1046.5-fqcy.smr'};...
     {'c1054.2-fqcy.smr'};{'c1307-fqcy.smr'};{'c1520.9-fqcy.smr'};...
     {'c2051_fqcy.smr'};{'c1176.3-evokfqcy1.smr'};{'c1277.6-fqcy2evok.smr'};...
     {'c1615.5-fqcyevok.smr'};{'c2098.4-fqcyEvok.smr'}];
 v_filesSpont = [{'c976.4-fqcy.smr'};{'c1036.2_fqcy.smr'};{'c1046.5-fqcy.smr'};...
     {'c1054.2-fqcy.smr'};{'c1307-fqcy.smr'};{'c1520.9-fqcy.smr'};...
     {'c2051_fqcy.smr'};{'c1176.3-spontfqcy.smr'};{'c1277.6-fqcySpont.smr'};...
     {'c1615.5-fqcyspont.smr'};{'c2098.4-fqcySpont.smr'}];

s_nbFiles = length(v_filesEvok);
 
m_corrsSpont = repmat(struct('number',0,'length',0,'corr_sp',0,'timelag_sp',0,...
    'corr_gam',0,'timelag_gam',0,'corr_gamLo',0,'timelag_gamLo',0,...
    'corr_gamMid',0,'timelag_gamMid',0,'corr_gamHi',0,'timelag_gamHi',0)...
    ,s_nbFiles*20,1);
m_corrsEvok = repmat(struct('number',0,'length',0,'corr_sp',0,'timelag_sp',0,...
    'corr_gam',0,'timelag_gam',0,'corr_gamLo',0,'timelag_gamLo',0,...
    'corr_gamMid',0,'timelag_gamMid',0,'corr_gamHi',0,'timelag_gamHi',0)...
    ,s_nbFiles*20,1);

s_lag = 0;
for s_file = 1:s_nbFiles
    st_info = openSpike2(v_filesEvok{s_file});
    [~,m_evok,~] = spont_evok_extract(st_info,false,true);
    st_info = openSpike2(v_filesSpont{s_file});
    [m_spont,~,st_info] = spont_evok_extract(st_info,true,false);
    s_SI = st_info(1).header.sampleinterval;
    s_Fe = st_info(1).header.Sampling;
    
%% Spontaneous upstates %%
    for s_which = 1:length(m_spont)
        v_Vm = m_spont{s_which,1};
        v_EEG = m_spont{s_which,2};
        v_diffVm = (v_Vm(2:end) - v_Vm(1:end-1))/360*1000;

%% select "good" upstates and filter %%
        b_valid = validity_test(v_Vm, st_info(1).header.sampleinterval);
        if b_valid
            [v_FiltSpVm,v_FiltGamLoVm,v_FiltGamMidVm, v_FiltGamHiVm,...
            v_FiltSpEEG,v_FiltGamLoEEG,v_FiltGamMidEEG,v_FiltGamHiEEG] =...
            filtre_upstate(v_Vm,v_EEG,st_info(1).header.Sampling);
            assert((length(v_FiltSpVm)==length(v_FiltGamLoVm)) && ...
                    (length(v_FiltSpVm)==length(v_FiltGamMidVm)) && ...
                    (length(v_FiltSpVm)==length(v_FiltGamHiVm)) && ...
                    (length(v_FiltSpVm)==length(v_FiltSpEEG)) && ...
                    (length(v_FiltSpEEG)==length(v_FiltGamLoEEG)) && ...
                    (length(v_FiltSpEEG)==length(v_FiltGamMidEEG)) && ...
                    (length(v_FiltSpEEG)==length(v_FiltGamHiEEG)));
            s_len = length(v_FiltSpVm);
            m_corrsSpont(s_which+s_lag).number = s_which;
            m_corrsSpont(s_which+s_lag).length = s_len*0.36;

%% compute correlations and timelag %%
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
            m_corrsSpont(s_which+s_lag).corr_sp = s_corrSp;
            m_corrsSpont(s_which+s_lag).timelag_sp = s_timeDiffSp;

            [s_corrGamLo,I_gamLo] = max(abs(v_xcorrGamLo));
            s_timeDiffGamLo = v_lagGamLo(I_gamLo)/s_Fe;
            m_corrsSpont(s_which+s_lag).corr_gamLo = s_corrGamLo;
            m_corrsSpont(s_which+s_lag).timelag_gamLo = s_timeDiffGamLo;

            [s_corrGamMid,I_gamMid] = max(abs(v_xcorrGamMid));
            s_timeDiffGamMid = v_lagGamMid(I_gamMid)/s_Fe;
            m_corrsSpont(s_which+s_lag).corr_gamMid = s_corrGamMid;
            m_corrsSpont(s_which+s_lag).timelag_gamMid = s_timeDiffGamMid;

            [s_corrGamHi,I_gamHi] = max(abs(v_xcorrGamHi));
            s_timeDiffGamHi = v_lagGamHi(I_gamHi)/s_Fe;
            m_corrsSpont(s_which+s_lag).corr_gamHi = s_corrGamHi;
            m_corrsSpont(s_which+s_lag).timelag_gamHi = s_timeDiffGamHi;

        else continue
        end
    end

%% Evoked upstates %%
    for s_which=1:length(m_evok)

        v_Vm = m_evok{s_which,1};
        v_EEG = m_evok{s_which,2};
        v_diffVm = (v_Vm(2:end) - v_Vm(1:end-1))/360*1000;

%% select "good" upstates and filter %%
        b_len = length(v_Vm)*s_SI/1000 > 300;
        b_pa = sum(v_diffVm >10) <0.5;
        if b_len && b_pa
            [v_FiltSpVm,v_FiltGamLoVm,v_FiltGamMidVm, v_FiltGamHiVm,...
            v_FiltSpEEG,v_FiltGamLoEEG,v_FiltGamMidEEG,v_FiltGamHiEEG] =...
            filtre_upstate(v_Vm,v_EEG,st_info(1).header.Sampling);
            assert((length(v_FiltSpVm)==length(v_FiltGamLoVm)) && ...
                    (length(v_FiltSpVm)==length(v_FiltGamMidVm)) && ...
                    (length(v_FiltSpVm)==length(v_FiltGamHiVm)) && ...
                    (length(v_FiltSpVm)==length(v_FiltSpEEG)) && ...
                    (length(v_FiltSpEEG)==length(v_FiltGamLoEEG)) && ...
                    (length(v_FiltSpEEG)==length(v_FiltGamMidEEG)) && ...
                    (length(v_FiltSpEEG)==length(v_FiltGamHiEEG)));
            s_len = length(v_FiltSpVm);
            m_corrsEvok(s_which+s_lag).number = s_which;
            m_corrsEvok(s_which+s_lag).length = s_len*0.36;

 %% Correlations and time lags %%
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
            lagDiff_sp = v_lagSp(I_sp);s_timeDiffSp = lagDiff_sp/s_Fe;
            m_corrsEvok(s_which+s_lag).corr_sp = s_corrSp;
            m_corrsEvok(s_which+s_lag).timelag_sp = s_timeDiffSp;

            [s_corrGamLo,I_gamLo] = max(abs(v_xcorrGamLo));
            lagDiff_gamLo = v_lagGamLo(I_gamLo);s_timeDiffGamLo = lagDiff_gamLo/s_Fe;
            m_corrsEvok(s_which+s_lag).corr_gamLo = s_corrGamLo;
            m_corrsEvok(s_which+s_lag).timelag_gamLo = s_timeDiffGamLo;

            [s_corrGamMid,I_gamMid] = max(abs(v_xcorrGamMid));
            lagDiff_gamMid = v_lagGamMid(I_gamMid);s_timeDiffGamMid = lagDiff_gamMid/s_Fe;
            m_corrsEvok(s_which+s_lag).corr_gamMid = s_corrGamMid;
            m_corrsEvok(s_which+s_lag).timelag_gamMid = s_timeDiffGamMid;

            [s_corrGamHi,I_gamHi] = max(abs(v_xcorrGamHi));
            lagDiff_gamHi = v_lagGamHi(I_gamHi);s_timeDiffGamHi = lagDiff_gamHi/s_Fe;
            m_corrsEvok(s_which+s_lag).corr_gamHi = s_corrGamHi;
            m_corrsEvok(s_which+s_lag).timelag_gamHi = s_timeDiffGamHi;

        else continue
        end
    end
    s_lag = s_lag + length(m_spont);
end