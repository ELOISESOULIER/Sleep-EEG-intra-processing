                                                                                                                                                                                                                                
v_filesEvok = [{'c976.4-fqcy.smr'};{'c1046.5-fqcy.smr'};...
     {'c1054.2-fqcy.smr'};{'c1307-fqcy.smr'};{'c1520.9-fqcy.smr'};...
     {'c2051_fqcy.smr'};{'c1615.5-fqcyevok.smr'};{'c2098.4-fqcyEvok.smr'}];
 v_filesSpont = [{'c976.4-fqcy.smr'};{'c1046.5-fqcy.smr'};...
     {'c1054.2-fqcy.smr'};{'c1307-fqcy.smr'};{'c1520.9-fqcy.smr'};...
     {'c2051_fqcy.smr'};{'c1615.5-fqcyspont.smr'};{'c2098.4-fqcySpont.smr'}];
 
assert(length(v_filesEvok) == length(v_filesSpont));
s_nbFiles = length(v_filesEvok);

params = getParams();

s_TFreso = params.s_TFreso;
s_nbBins = params.s_nbBins;
s_freqMin = params.s_freqMin;
s_sampInt =360;

m_TotAbsVmSp = [];
m_TotAbsVmEv = [];
m_TotDiffAbsVm = [];
m_TotAbsEEGSp = [];
m_TotAbsEEGEv = [];
m_TotDiffAbsEEG = [];

m_TotAbsVmSpNorm = [];
m_TotAbsVmEvNorm = [];
m_TotDiffAbsVmNorm = [];
m_TotAbsEEGSpNorm = [];
m_TotAbsEEGEvNorm = [];
m_TotDiffAbsEEGNorm = [];

m_TotRelVmSp = [];
m_TotRelVmEv = [];
m_TotDiffRelVm = [];
m_TotRelEEGSp = [];
m_TotRelEEGEv = [];
m_TotDiffRelEEG = [];

m_TotRelVmSpNorm = [];
m_TotRelVmEvNorm = [];
m_TotDiffRelVmNorm = [];
m_TotRelEEGSpNorm = [];
m_TotRelEEGEvNorm = [];
m_TotDiffRelEEGNorm = [];

v_TotRawVmSp = [];
v_TotSpVmSp = [];
v_TotGamLoVmSp = [];
v_TotGamMidVmSp = [];
v_TotGamHiVmSp = [];

v_TotRawEEGSp = [];
v_TotSpEEGSp = [];
v_TotGamLoEEGSp = [];
v_TotGamMidEEGSp = [];
v_TotGamHiEEGSp = [];

v_TotRawVmEv = [];
v_TotSpVmEv = [];
v_TotGamLoVmEv = [];
v_TotGamMidVmEv = [];
v_TotGamHiVmEv = [];

v_TotRawEEGEv = [];
v_TotSpEEGEv = [];
v_TotGamLoEEGEv = [];
v_TotGamMidEEGEv = [];
v_TotGamHiEEGEv = [];

s_countSpont = 0;
s_countEvok = 0;    

 for index = 1:s_nbFiles
    
    %check whether spontaneous and evoked upstates are on the same file 
    b_divided = v_filesEvok{index}~=v_filesSpont{index}(1:length(v_filesEvok{index}));

    if sum(b_divided)>1
        st_infoSpont = openSpike2(v_filesSpont{index});
        st_infoEvok = openSpike2(v_filesEvok{index});
    else
        st_info = openSpike2(v_filesSpont{index});
    end

    v_FreqAxis = linspace(100,s_freqMin,s_TFreso);

     %% Time-Frequency maps of all kinds, relatively averaged or not, normalized or not
     
     if sum(b_divided)>1 %compute everything separately for evoked and spontaneous conditions
         
        [m_spont,~,st_infoSpont] = spont_evok_extract(st_infoSpont,true,false,true);
        [~,m_evok,st_infoEvok] = spont_evok_extract(st_infoEvok,false,true,true); 
        
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
        
%% Compute map of differences
        m_DiffAbsVm = m_TFAbsVmSpont - m_TFAbsVmEvok;
        m_DiffAbsEEG = m_TFAbsEEGSpont - m_TFAbsEEGEvok;
        m_DiffRelVm = m_TFRelVmSpont - m_TFRelVmEvok;
        m_DiffRelEEG = m_TFRelEEGSpont - m_TFRelEEGEvok;
        
        m_DiffAbsVmNorm = m_TFAbsVmSpontNorm - m_TFAbsVmEvokNorm;
        m_DiffAbsEEGNorm = m_TFAbsEEGSpontNorm - m_TFAbsEEGEvokNorm;
        m_DiffRelVmNorm = m_TFRelVmSpontNorm - m_TFRelVmEvokNorm;
        m_DiffRelEEGNorm = m_TFRelEEGSpontNorm - m_TFRelEEGEvokNorm;
        
        
     else %compute them together
        [m_spont,m_evok,st_info] = spont_evok_extract(st_info,true,true,true);

        [m_TFAbsVmSpont,m_TFAbsEEGSpont,m_TFAbsVmEvok,m_TFAbsEEGEvok,v_TimeAxis]=...
             TFabs1sec_spontevok(st_info,true,true,s_TFreso,false);
        [m_TFRelVmSpont,m_TFRelEEGSpont,m_TFRelVmEvok,m_TFRelEEGEvok]=...
            TFrelatif_spontevok(st_info,true,true,s_nbBins,s_TFreso,false);

        [m_TFAbsVmSpontNorm,m_TFAbsEEGSpontNorm,m_TFAbsVmEvokNorm,...
            m_TFAbsEEGEvokNorm] = TFabs1sec_spontevok(st_info);
        [m_TFRelVmSpontNorm,m_TFRelEEGSpontNorm,m_TFRelVmEvokNorm,m_TFRelEEGEvokNorm]=...
            TFrelatif_spontevok(st_info);   
        
        % Compute map of differences
        m_DiffAbsVm = m_TFAbsVmSpont - m_TFAbsVmEvok;
        m_DiffAbsEEG = m_TFAbsEEGSpont - m_TFAbsEEGEvok;
        m_DiffRelVm = m_TFRelVmSpont - m_TFRelVmEvok;
        m_DiffRelEEG = m_TFRelEEGSpont - m_TFRelEEGEvok;
        
        m_DiffAbsVmNorm = m_TFAbsVmSpontNorm - m_TFAbsVmEvokNorm;
        m_DiffAbsEEGNorm = m_TFAbsEEGSpontNorm - m_TFAbsEEGEvokNorm;
        m_DiffRelVmNorm = m_TFRelVmSpontNorm - m_TFRelVmEvokNorm;
        m_DiffRelEEGNorm = m_TFRelEEGSpontNorm - m_TFRelEEGEvokNorm;

     end
    
     
     
     figure(6)
     subplot(2,2,1)
     f_ImageArray(m_DiffAbsVm,v_TimeAxis,v_FreqAxis,'colormap','jet');
     xlabel('Temps (s)');ylabel('Frequence (Hz)');title('Absolu-Intra')
     subplot(2,2,2)
     f_ImageArray(m_DiffAbsVmNorm,v_TimeAxis,v_FreqAxis,'colormap','jet');
     xlabel('Temps (s)');title('Absolu-Intra-Norm');
     subplot(2,2,3)
     f_ImageArray(m_DiffRelVm,1:1:s_nbBins,v_FreqAxis,'colormap','jet');
     xlabel('Temps (%upstate)');ylabel('Fréquence (Hz)');title('Relatif-Intra');
     subplot(2,2,4)
     f_ImageArray(m_DiffRelVmNorm,1:1:s_nbBins,v_FreqAxis,'colormap','jet');
     xlabel('Temps (%upstate)');ylabel('Fréquence (Hz)'); title('Relatif-Intra-Norm');
     
     figure(7)
     subplot(2,2,1)
     f_ImageArray(m_DiffAbsEEG,v_TimeAxis,v_FreqAxis,'colormap','jet');
     xlabel('Temps (s)');ylabel('Frequence (Hz)');title('Absolu-EEG')
     subplot(2,2,2)
     f_ImageArray(m_DiffAbsEEGNorm,v_TimeAxis,v_FreqAxis,'colormap','jet');
     xlabel('Temps (s)');title('Absolu-EEG-Norm');
     subplot(2,2,3)
     f_ImageArray(m_DiffRelEEG,1:1:s_nbBins,v_FreqAxis,'colormap','jet');
     xlabel('Temps (%upstate)');ylabel('Fréquence (Hz)');title('Relatif-EEG');
     subplot(2,2,4)
     f_ImageArray(m_DiffRelEEGNorm,1:1:s_nbBins,v_FreqAxis,'colormap','jet');
     xlabel('Temps (%upstate)');ylabel('Fréquence (Hz)'); title('Relatif-EEG-Norm');
     
    %% Compute mean time frequency maps
    
    % Abs & no normalization
    if isempty(m_TotAbsVmSp)
        m_TotAbsVmSp = m_TFAbsVmSpont;
    else
        m_TotAbsVmSp = m_TotAbsVmSp + m_TFAbsVmSpont;
    end
    if isempty(m_TotAbsVmEv)
        m_TotAbsVmEv = m_TFAbsVmEvok;
    else
        m_TotAbsVmEv = m_TotAbsVmEv + m_TFAbsVmEvok;
    end
    if isempty(m_TotAbsEEGSp)
        m_TotAbsEEGSp = m_TFAbsEEGSpont;
    else
        m_TotAbsEEGSp = m_TotAbsEEGSp + m_TFAbsEEGSpont;
    end
     if isempty(m_TotAbsEEGEv)
        m_TotAbsEEGEv = m_TFAbsEEGEvok;
    else
        m_TotAbsEEGEv = m_TotAbsEEGEv + m_TFAbsEEGEvok;
     end
     if isempty(m_TotDiffAbsVm)
         m_TotDiffAbsVm = m_DiffAbsVm;
     else m_TotDiffAbsVm = m_TotDiffAbsVm + m_DiffAbsVm;
     end
     if isempty(m_TotDiffAbsEEG)
         m_TotDiffAbsEEG = m_DiffAbsEEG;
     else m_TotDiffAbsEEG = m_TotDiffAbsEEG + m_DiffAbsEEG;
     end
     
     %Abs & normalization
     if isempty(m_TotAbsVmSpNorm)
        m_TotAbsVmSpNorm = m_TFAbsVmSpontNorm;
    else
        m_TotAbsVmSpNorm = m_TotAbsVmSpNorm + m_TFAbsVmSpontNorm;
     end
    if isempty(m_TotAbsVmEvNorm)
        m_TotAbsVmEvNorm = m_TFAbsVmEvokNorm;
    else
        m_TotAbsVmEvNorm = m_TotAbsVmEvNorm + m_TFAbsVmEvokNorm;
    end
    if isempty(m_TotAbsEEGSpNorm)
        m_TotAbsEEGSpNorm = m_TFAbsEEGSpontNorm;
    else
        m_TotAbsEEGSpNorm = m_TotAbsEEGSpNorm + m_TFAbsEEGSpontNorm;
    end
     if isempty(m_TotAbsEEGEvNorm)
        m_TotAbsEEGEvNorm = m_TFAbsEEGEvokNorm;
    else
        m_TotAbsEEGEvNorm = m_TotAbsEEGEvNorm + m_TFAbsEEGEvokNorm;
     end
     if isempty(m_TotDiffAbsVmNorm)
         m_TotDiffAbsVmNorm = m_DiffAbsVmNorm;
     else m_TotDiffAbsVmNorm = m_TotDiffAbsVmNorm + m_DiffAbsVmNorm;
     end
     if isempty(m_TotDiffAbsEEGNorm)
         m_TotDiffAbsEEGNorm = m_DiffAbsEEGNorm;
     else m_TotDiffAbsEEGNorm = m_TotDiffAbsEEGNorm + m_DiffAbsEEGNorm;
     end
    
      % Rel & no normalization
    if isempty(m_TotRelVmSp)
        m_TotRelVmSp = m_TFRelVmSpont;
    else
        m_TotRelVmSp = m_TotRelVmSp + m_TFRelVmSpont;
    end
    if isempty(m_TotRelVmEv)
        m_TotRelVmEv = m_TFRelVmEvok;
    else
        m_TotRelVmEv = m_TotRelVmEv + m_TFRelVmEvok;
    end
    if isempty(m_TotRelEEGSp)
        m_TotRelEEGSp = m_TFRelEEGSpont;
    else
        m_TotRelEEGSp = m_TotRelEEGSp + m_TFRelEEGSpont;
    end
     if isempty(m_TotRelEEGEv)
        m_TotRelEEGEv = m_TFRelEEGEvok;
    else
        m_TotRelEEGEv = m_TotRelEEGEv + m_TFRelEEGEvok;
     end
     
     %Rel & normalization
     if isempty(m_TotRelVmSpNorm)
        m_TotRelVmSpNorm = m_TFRelVmSpontNorm;
    else
        m_TotRelVmSpNorm = m_TotRelVmSpNorm + m_TFRelVmSpontNorm;
     end
    if isempty(m_TotRelVmEvNorm)
        m_TotRelVmEvNorm = m_TFRelVmEvokNorm;
    else
        m_TotRelVmEvNorm = m_TotRelVmEvNorm + m_TFRelVmEvokNorm;
    end
    if isempty(m_TotRelEEGSpNorm)
        m_TotRelEEGSpNorm = m_TFRelEEGSpontNorm;
    else
        m_TotRelEEGSpNorm = m_TotRelEEGSpNorm + m_TFRelEEGSpontNorm;
    end
     if isempty(m_TotRelEEGEvNorm)
        m_TotRelEEGEvNorm = m_TFRelEEGEvokNorm;
    else
        m_TotRelEEGEvNorm = m_TotRelEEGEvNorm + m_TFRelEEGEvokNorm;
     end
    
%       
    %% normalization factor for EEG
        s_meanEEG = 0;
        s_stdEEG = 0; 
        s_nbEEG = 0;
        for s_which = 1:length(m_spont)
            b_valid = validity_test(m_spont{s_which,1},st_info(1).header.sampleinterval);
            if b_valid
                s_meanEEG = s_meanEEG + mean(m_spont{s_which,2});
                s_stdEEG = s_stdEEG + std(m_spont{s_which,2});
                s_nbEEG = s_nbEEG + 1;
            end
        end
        s_meanEEG = s_meanEEG/s_nbEEG;
        s_stdEEG = s_stdEEG/s_nbEEG;

    %% Raw and filtered signal    
%     for s_which=1:length(m_spont)
%         [b_valid,m_spont{s_which,1}] = validity_test(m_spont{s_which,1},...
%             st_info(1).header.sampleinterval);
%         if b_valid
%             s_countSpont = s_countSpont + 1;
%             m_spont{s_which,1} = m_spont{s_which,1} -70 - m_spont{s_which,1}(1);
%             m_spont{s_which,2} = (m_spont{s_which,2}-s_meanEEG)/s_stdEEG^2;
%             [v_FiltSpVm,v_FiltGamLoVm,v_FiltGamMidVm, v_FiltGamHiVm,...
%                 v_FiltSpEEG,v_FiltGamLoEEG,v_FiltGamMidEEG,v_FiltGamHiEEG] =...
%                 filtre_upstate(m_spont{s_which,1},m_spont{s_which,2},st_info(1).header.Sampling);
%             
%             if isempty(v_TotRawVmSp)
%             v_TotRawVmSp = m_spont{s_which,1};
%             else v_TotRawVmSp = v_TotRawVmSp + m_spont{s_which,1};
%             end
%             if isempty(v_TotSpVmSp)
%                 v_TotSpVmSp = v_FiltSpVm;
%             else v_TotSpVmSp = v_TotSpVmSp + v_FiltSpVm;
%             end
%             if isempty(v_TotGamLoVmSp)
%                 v_TotGamLoVmSp = v_FiltGamLoVm;
%             else v_TotGamLoVmSp = v_TotGamLoVmSp + v_FiltGamLoVm;
%             end
%             if isempty(v_TotGamMidVmSp)
%                 v_TotGamMidVmSp = v_FiltGamMidVm;
%             else v_TotGamMidVmSp = v_TotGamMidVmSp + v_FiltGamMidVm;
%             end
%             if isempty(v_TotGamHiVmSp)
%                 v_TotGamHiVmSp = v_FiltGamHiVm;
%             else v_TotGamHiVmSp = v_TotGamHiVmSp + v_FiltGamHiVm;
%             end
%             
%             if isempty(v_TotRawEEGSp)
%                 v_TotRawEEGSp = m_spont{s_which,2};
%             else v_TotRawEEGSp = v_TotRawEEGSp + m_spont{s_which,2};
%             end
%             if isempty(v_TotSpEEGSp)
%                 v_TotSpEEGSp = v_FiltSpEEG;                     
%             else v_TotSpEEGSp = v_TotSpEEGSp + v_FiltSpEEG;
%             end
%             if isempty(v_TotGamLoEEGSp)
%                 v_TotGamLoEEGSp = v_FiltGamLoEEG;
%             else v_TotGamLoEEGSp = v_TotGamLoEEGSp + v_FiltGamLoEEG;
%             end
%             if isempty(v_TotGamMidEEGSp)
%                 v_TotGamMidEEGSp = v_FiltGamMidEEG;
%             else v_TotGamMidEEGSp = v_TotGamMidEEGSp + v_FiltGamMidEEG;
%             end
%             if isempty(v_TotGamHiEEGSp)
%                 v_TotGamHiEEGSp = v_FiltGamHiEEG;
%             else v_TotGamHiEEGSp = v_TotGamHiEEGSp + v_FiltGamHiEEG;
%             end
%             
%         end
%     end
% 
%     for s_which=1:length(m_evok)
%         [b_valid,m_evok{s_which,1}] = validity_test(m_evok{s_which,1},...
%             st_info(1).header.sampleinterval,true);
%         if b_valid
%             s_countEvok = s_countEvok + 1;
%             m_evok{s_which,1} = m_evok{s_which,1}-70 -m_evok{s_which,1}(1);
%             m_evok{s_which,2} = (m_evok{s_which,2}-s_meanEEG)/s_stdEEG^2;
%             [v_FiltSpVm,v_FiltGamLoVm,v_FiltGamMidVm, v_FiltGamHiVm,...
%                 v_FiltSpEEG,v_FiltGamLoEEG,v_FiltGamMidEEG,v_FiltGamHiEEG] =...
%                 filtre_upstate(m_evok{s_which,1},m_evok{s_which,2},st_info(1).header.Sampling);
%        
%             if isempty(v_TotRawVmEv)
%                 v_TotRawVmEv = m_evok{s_which,1};
%             else v_TotRawVmEv = v_TotRawVmEv + m_evok{s_which,1};
%             end
%             if isempty(v_TotSpVmEv)
%                 v_TotSpVmEv = v_FiltSpVm;
%             else v_TotSpVmEv = v_TotSpVmEv + v_FiltSpVm;
%             end
%             if isempty(v_TotGamLoVmEv)
%                 v_TotGamLoVmEv = v_FiltGamLoVm;
%             else v_TotGamLoVmEv = v_TotGamLoVmEv + v_FiltGamLoVm;
%             end
%             if isempty(v_TotGamMidVmEv)
%                 v_TotGamMidVmEv = v_FiltGamMidVm;
%             else v_TotGamMidVmEv = v_TotGamMidVmEv + v_FiltGamMidVm;
%             end
%             if isempty(v_TotGamHiVmEv)
%                 v_TotGamHiVmEv = v_FiltGamHiVm;
%             else v_TotGamHiVmEv = v_TotGamHiVmEv + v_FiltGamHiVm;
%             end
%             
%             if isempty(v_TotRawEEGEv)
%                 v_TotRawEEGEv = m_evok{s_which,2};
%             else v_TotRawEEGEv = v_TotRawEEGEv + m_evok{s_which,2};
%             end
%             if isempty(v_TotSpEEGEv)
%                 v_TotSpEEGEv = v_FiltSpEEG;
%             else v_TotSpEEGEv = v_TotSpEEGEv + v_FiltSpEEG;
%             end
%             if isempty(v_TotGamLoEEGEv)
%                 v_TotGamLoEEGEv = v_FiltGamLoEEG;
%             else v_TotGamLoEEGEv = v_TotGamLoEEGEv + v_FiltGamLoEEG;
%             end
%             if isempty(v_TotGamMidEEGEv)
%                 v_TotGamMidEEGEv = v_FiltGamMidEEG;
%             else v_TotGamMidEEGEv = v_TotGamMidEEGEv + v_FiltGamMidEEG;
%             end
%             if isempty(v_TotGamHiEEGEv)
%                 v_TotGamHiEEGEv = v_FiltGamHiEEG;
%             else v_TotGamHiEEGEv = v_TotGamHiEEGEv + v_FiltGamHiEEG;
%             end
%             
%         end
%     end
 end

%% Average all

m_TotAbsVmSp = m_TotAbsVmSp/s_nbFiles;
m_TotAbsVmEv = m_TotAbsVmEv/s_nbFiles;
m_TotDiffAbsVm = m_TotDiffAbsVm/s_nbFiles;
m_TotAbsEEGSp = m_TotAbsEEGSp/s_nbFiles;
m_TotAbsEEGEv = m_TotAbsEEGEv/s_nbFiles;
m_TotDiffAbsEEG = m_TotDiffAbsEEG/s_nbFiles;

m_TotAbsVmSpNorm = m_TotAbsVmSpNorm/s_nbFiles;
m_TotAbsVmEvNorm = m_TotAbsVmEvNorm/s_nbFiles;
m_TotDiffAbsVmNorm = m_TotDiffAbsVmNorm/s_nbFiles;
m_TotAbsEEGSpNorm = m_TotAbsEEGSpNorm/s_nbFiles;
m_TotAbsEEGEvNorm = m_TotAbsEEGEvNorm/s_nbFiles;
m_TotDiffAbsEEGNorm = m_TotDiffAbsEEGNorm/s_nbFiles;

m_TotRelVmSp = m_TotRelVmSp/s_nbFiles;
m_TotRelVmEv = m_TotRelVmEv/s_nbFiles;
m_TotDiffRelVm = m_TotDiffRelVm/s_nbFiles;
m_TotRelEEGSp = m_TotRelEEGSp/s_nbFiles;
m_TotRelEEGEv = m_TotRelEEGEv/s_nbFiles;
m_TotDiffRelEEG = m_TotDiffRelEEG/s_nbFiles;

m_TotRelVmSpNorm = m_TotRelVmSpNorm/s_nbFiles;
m_TotRelVmEvNorm = m_TotRelVmEvNorm/s_nbFiles;
m_TotDiffRelVmNorm = m_TotDiffRelVmNorm/s_nbFiles;
m_TotRelEEGSpNorm = m_TotRelEEGSpNorm/s_nbFiles;
m_TotRelEEGEvNorm = m_TotRelEEGEvNorm/s_nbFiles;
m_TotDiffRelEEGNorm = m_TotDiffRelEEGNorm/s_nbFiles;

v_TotRawVmSp = v_TotRawVmSp/s_countSpont;
v_TotSpVmSp = v_TotSpVmSp/s_countSpont;
v_TotGamLoVmSp = v_TotGamLoVmSp/s_countSpont;
v_TotGamMidVmSp = v_TotGamMidVmSp/s_countSpont;
v_TotGamHiVmSp = v_TotGamHiVmSp/s_countSpont;

v_TotRawEEGSp = v_TotRawEEGSp/s_countSpont;
v_TotSpEEGSp = v_TotSpEEGSp/s_countSpont;
v_TotGamLoEEGSp = v_TotGamLoEEGSp/s_countSpont;
v_TotGamMidEEGSp = v_TotGamMidEEGSp/s_countSpont;
v_TotGamHiEEGSp = v_TotGamHiEEGSp/s_countSpont;

v_TotRawVmEv = v_TotRawVmEv/s_countEvok;
v_TotSpVmEv = v_TotSpVmEv/s_countEvok;
v_TotGamLoVmEv = v_TotGamLoVmEv/s_countEvok;
v_TotGamMidVmEv = v_TotGamMidVmEv/s_countEvok;
v_TotGamHiVmEv = v_TotGamHiVmEv/s_countEvok;

v_TotRawEEGEv = v_TotRawEEGEv/s_countEvok;
v_TotSpEEGEv = v_TotSpEEGEv/s_countEvok;
v_TotGamLoEEGEv = v_TotGamLoEEGEv/s_countEvok;
v_TotGamMidEEGEv = v_TotGamMidEEGEv/s_countEvok;
v_TotGamHiEEGEv = v_TotGamHiEEGEv/s_countEvok;

s_cut = floor(400/s_sampInt*1000);
v_TimeAxis = v_TimeAxis(1:end-s_cut+1);

m_TotAbsVmSpNorm = m_TotAbsVmSpNorm(:,s_cut:end);
m_TotAbsVmEvNorm = m_TotAbsVmEvNorm(:,s_cut:end);
m_TotDiffAbsVmNorm = m_TotDiffAbsVmNorm(:,s_cut:end);

m_TotAbsEEGSpNorm = m_TotAbsEEGSpNorm(:,s_cut:end);
m_TotAbsEEGEvNorm = m_TotAbsEEGEvNorm(:,s_cut:end);
m_TotDiffAbsEEGNorm = m_TotDiffAbsEEGNorm(:,s_cut:end);

m_TotAbsVmSp = m_TotAbsVmSp(:,s_cut:end);
m_TotAbsVmEv = m_TotAbsVmEv(:,s_cut:end);
m_TotDiffAbsVm = m_TotDiffAbsVm(:,s_cut:end);

m_TotAbsEEGSp = m_TotAbsEEGSp(:,s_cut:end);
m_TotAbsEEGEv = m_TotAbsEEGEv(:,s_cut:end);
m_TotDiffAbsEEG = m_TotDiffAbsEEG(:,s_cut:end);

%% Visualize all

figure(20)
subplot(2,2,1)
    v_clim = [prctile(m_TotAbsVmSpNorm(:),3),prctile(m_TotAbsVmSpNorm(:),100)];
    f_ImageArray(m_TotAbsVmSpNorm,v_TimeAxis,v_FreqAxisTot,'colormap','jet','limits',v_clim);
    title('TF absolu - Intra - Spontané');xlabel('Temps (s)');ylabel('Frequence (Hz)')
subplot(2,2,2)
  v_clim = [prctile(m_TotAbsVmEvNorm(:),3),prctile(m_TotAbsVmEvNorm(:),100)];
    f_ImageArray(m_TotAbsVmEvNorm,v_TimeAxis,v_FreqAxisTot,'colormap','jet','limits',v_clim);
    title('TF absolu - Intra - Évoqué');xlabel('Temps (s)');ylabel('Frequence (Hz)')
subplot(2,2,3)
     v_clim = [prctile(m_TotRelVmSpNorm(:),1),prctile(m_TotRelVmSpNorm(:),99)];
     f_ImageArray(m_TotRelVmSpNorm,(1:1:s_nbBins)*1/s_nbBins,v_FreqAxisTot,'colormap','jet','limits',v_clim);
     title('TF relatif - Intra - Spontané');xlabel('Temps (%upstate)');ylabel('Frequence (Hz)')
subplot(2,2,4)
     v_clim = [prctile(m_TotRelVmEvNorm(:),1),prctile(m_TotRelVmEvNorm(:),99)];
     f_ImageArray(m_TotRelVmEvNorm,(1:1:s_nbBins)*1/s_nbBins,v_FreqAxisTot,'colormap','jet','limits',v_clim);
     title('TF relatif - Intra - Évoqué');xlabel('Temps (%upstate)');ylabel('Frequence (Hz)')
    
figure(21)
subplot(2,2,1)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
    v_clim = [prctile(m_TotAbsEEGSpNorm(:),5),prctile(m_TotAbsEEGSpNorm(:),100)];
    f_ImageArray(m_TotAbsEEGSpNorm,v_TimeAxis,v_FreqAxisTot,'colormap','jet','limits',v_clim);
    title('TF absolu - EEG - Spontané');xlabel('Temps (s)');ylabel('Frequence (Hz)')                                                                     
subplot(2,2,2)
    v_clim = [prctile(m_TotAbsEEGEvNorm(:),5),prctile(m_TotAbsEEGEvNorm(:),100)];
    f_ImageArray(m_TotAbsEEGEvNorm,v_TimeAxis,v_FreqAxisTot,'colormap','jet','limits',v_clim);
    title('TF absolu - EEG - Évoqué');xlabel('Temps (s)');ylabel('Frequence (Hz)')
subplot(2,2,3)
    v_clim = [prctile(m_TotRelEEGSpNorm(:),1),prctile(m_TotRelEEGSpNorm(:),100)];
    f_ImageArray(m_TotRelEEGSpNorm,(1:1:s_nbBins)*1/s_nbBins,v_FreqAxisTot,'colormap','jet','limits',v_clim);
    title('TF relatif - EEG - Spontané');xlabel('Temps (%upstate)');ylabel('Frequence (Hz)')
subplot(2,2,4)
    v_clim = [prctile(m_TotRelEEGEvNorm(:),1),prctile(m_TotRelEEGEvNorm(:),100)];
    f_ImageArray(m_TotRelEEGEvNorm,(1:1:s_nbBins)*1/s_nbBins,v_FreqAxisTot,'colormap','jet','limits',v_clim);
    title('TF relatif - EEG - Évoqué');xlabel('Temps (%upstate)');ylabel('Frequence (Hz)')

figure(31)
subplot(2,2,1)
    v_clim = [prctile(m_TotAbsVmSp(:),1),prctile(m_TotAbsVmSp(:),99)];
    f_ImageArray(m_TotAbsVmSp,v_TimeAxis,v_FreqAxis,'colormap','jet','limits',v_clim);
    title('TF absolu - Intra - Spontané');xlabel('Temps (s)');ylabel('Frequence (Hz)')
subplot(2,2,2)
    v_clim = [prctile(m_TotAbsVmEv(:),1),prctile(m_TotAbsVmEv(:),99)];
    f_ImageArray(m_TotAbsVmEv,v_TimeAxis,v_FreqAxis,'colormap','jet','limits',v_clim);
    title('TF absolu - Intra - Évoqué');xlabel('Temps (s)');ylabel('Frequence (Hz)')
subplot(2,2,3)
     v_clim = [prctile(m_TotRelVmSp(:),1),prctile(m_TotRelVmSp(:),99)];
     f_ImageArray(m_TotRelVmSp,(1:1:s_nbBins)*1/s_nbBins,v_FreqAxis,'colormap','jet','limits',v_clim);
     title('TF relatif - Intra - Spontané');xlabel('Temps (%upstate)');ylabel('Frequence (Hz)')
subplot(2,2,4)
     v_clim = [prctile(m_TotRelVmEv(:),1),prctile(m_TotRelVmEv(:),99)];
     f_ImageArray(m_TotRelVmEv,(1:1:s_nbBins)*1/s_nbBins,v_FreqAxis,'colormap','jet','limits',v_clim);
     title('TF relatif - Intra - Évoqué');xlabel('Temps (%upstate)');ylabel('Frequence (Hz)') 
    
figure(32)
subplot(2,2,1)
    v_clim = [prctile(m_TotAbsEEGSp(:),5),prctile(m_TotAbsEEGSp(:),100)];
    f_ImageArray(m_TotAbsEEGSp,v_TimeAxis,v_FreqAxis,'colormap','jet','limits',v_clim);
    title('TF absolu - EEG - Spontané');xlabel('Temps (s)');ylabel('Frequence (Hz)')
subplot(2,2,2)
    v_clim = [prctile(m_TotAbsEEGEv(:),5),prctile(m_TotAbsEEGEv(:),100)];
    f_ImageArray(m_TotAbsEEGEv,v_TimeAxis,v_FreqAxis,'colormap','jet','limits',v_clim);
    title('TF absolu - EEG - Évoqué');xlabel('Temps (s)');ylabel('Frequence (Hz)')
subplot(2,2,3)
    v_clim = [prctile(m_TotRelEEGSp(:),5),prctile(m_TotRelEEGSp(:),99)];
    f_ImageArray(m_TotRelEEGSp,(1:1:s_nbBins)*1/s_nbBins,v_FreqAxis,'colormap','jet','limits',v_clim);
    title('TF relatif - EEG - Spontané');xlabel('Temps (%upstate)');ylabel('Frequence (Hz)')
subplot(2,2,4)
    v_clim = [prctile(m_TotRelEEGEv(:),5),prctile(m_TotRelEEGEv(:),99)];
    f_ImageArray(m_TotRelEEGEv,(1:1:s_nbBins)*1/s_nbBins,v_FreqAxis,'colormap','jet','limits',v_clim);
    title('TF relatif - EEG - Évoqué');xlabel('Temps (%upstate)');ylabel('Frequence (Hz)')
 
figure(18)
    subplot(2,1,1)
     v_clim = [prctile(m_TotDiffAbsVm(:),1),prctile(m_TotDiffAbsVm(:),95)];
     f_ImageArray(m_TotDiffAbsVm,v_TimeAxis,v_FreqAxis,'colormap','jet','limits',v_clim);
     xlabel('Temps (s)');ylabel('Frequence (Hz)');title('Absolu-Intra')
     subplot(2,1,2)
     v_clim = [prctile(m_TotDiffAbsVmNorm(:),1),prctile(m_TotDiffAbsVmNorm(:),95)];
     f_ImageArray(m_TotDiffAbsVmNorm,v_TimeAxis,v_FreqAxis,'colormap','jet','limits',v_clim);
     xlabel('Temps (s)');title('Absolu-Intra-Norm');
%      subplot(2,2,3)
%      f_ImageArray(m_TotDiffRelVm,1:1:s_nbBins,v_FreqAxis,'colormap','jet');
%      xlabel('Temps (%upstate)');ylabel('Fréquence (Hz)');title('Relatif-Intra');
%      subplot(2,2,4)
%      f_ImageArray(m_TotDiffRelVmNorm,1:1:s_nbBins,v_FreqAxis,'colormap','jet');
%      xlabel('Temps (%upstate)');ylabel('Fréquence (Hz)'); title('Relatif-Intra-Norm');

figure(19)
    subplot(2,1,1)
     v_clim = [prctile(m_TotDiffAbsEEG(:),1),prctile(m_TotDiffAbsEEG(:),95)];
     f_ImageArray(m_TotDiffAbsEEG,v_TimeAxis,v_FreqAxis,'colormap','jet','limits',v_clim);
     xlabel('Temps (s)');ylabel('Frequence (Hz)');title('Absolu-EEG')
     subplot(2,1,2)
     v_clim = [prctile(m_TotDiffAbsEEGNorm(:),1),prctile(m_TotDiffAbsEEGNorm(:),95)];
     f_ImageArray(m_TotDiffAbsEEGNorm,v_TimeAxis,v_FreqAxis,'colormap','jet','limits',v_clim);
     xlabel('Temps (s)');title('Absolu-EEG-Norm');
     subplot(2,2,3)
     f_ImageArray(m_TotDiffRelEEG,1:1:s_nbBins,v_FreqAxis,'colormap','jet');
     xlabel('Temps (%upstate)');ylabel('Fréquence (Hz)');title('Relatif-Intra');
     subplot(2,2,4)
     f_ImageArray(m_TotDiffRelEEGNorm,1:1:s_nbBins,v_FreqAxis,'colormap','jet');
     xlabel('Temps (%upstate)');ylabel('Fréquence (Hz)'); title('Relatif-Intra-Norm');
     
   
%    figure(1)
%         time = 0:st_info(1).header.sampleinterval/1000000:...
%             st_info(1).header.sampleinterval/1000000*(length(m_spont{s_which,1})-1);
%         subplot(5,1,1);
%             plot(time, v_TotRawVmSp);hold on;
%             plot(time, v_TotRawVmEv);hold off;
%             legend('Spontané','Évoqué');
%             title('Intra'); ylabel('Potentiel (mV)')
%         subplot(5,1,2);
%             plot(time, v_TotSpVmSp(1:end-1));hold on;
%             plot(time,v_TotSpVmEv(1:end-1));hold off;
%             legend('Spontané','Évoqué');
%             title('Intra - Spindle');ylabel('Potentiel(mV)')
%         subplot(5,1,3)
%             plot(time,v_TotGamLoVmSp(1:end-1));hold on;
%             plot(time,v_TotGamLoVmEv(1:end-1));hold off;
%             legend('Spontané','Évoqué');
%             title('Intra - Gamma 30-60');ylabel('Potentiel(mV)')
%         subplot(5,1,4)
%             plot(time,v_TotGamMidVmSp(1:end-1));hold on;
%             plot(time,v_TotGamMidVmEv(1:end-1));hold off;
%             legend('Spontané','Évoqué');
%             title('Intra - Gamma 60-80');ylabel('Potentiel(mV)')
%         subplot(5,1,5)
%             plot(time,v_TotGamHiVmSp(1:end-1));hold on;
%             plot(time,v_TotGamHiVmEv(1:end-1));hold off;
%             legend('Spontané','Évoqué');
%             title('Intra - Gamma 80-100');ylabel('Potentiel(mV)')
%             
%        figure(2)     
%        subplot(5,1,1);
%             plot(time, v_TotRawEEGSp);hold on;
%             plot(time,v_TotRawEEGEv);hold off;
%             legend('Spontané','Évoqué');
%             title('EEG');ylabel('Potentiel (uV)')
%        subplot(5,1,2);
%             plot(time, v_TotSpEEGSp(1:end-1));hold on;
%             plot(time,v_TotSpEEGEv(1:end-1));hold off;
%             legend('Spontané','Évoqué');
%             title('EEG - Spindle');ylabel('Potentiel(uV)')
%        subplot(5,1,3);
%             plot(time,v_TotGamLoEEGSp(1:end-1));hold on;
%             plot(time,v_TotGamLoEEGEv(1:end-1));hold off;
%             legend('Spontané','Évoqué');
%             title('EEG - Gamma 30-60');ylabel('Potentiel (uV)')     
%        subplot(5,1,4);
%             plot(time,v_TotGamMidEEGSp(1:end-1));hold on;
%             plot(time,v_TotGamMidEEGEv(1:end-1));hold off;
%             legend('Spontané','Évoqué');
%             title('EEG - Gamma 60-80');ylabel('Potentiel(uV)')
%        subplot(5,1,5);
%             plot(time,v_TotGamHiEEGSp(1:end-1));hold on;
%             plot(time,v_TotGamHiEEGEv(1:end-1));hold on;
%             legend('Spontané','Évoqué');
%             title('EEG - Gamma 80-100');ylabel('Potentiel (uV)')
%   
