%computes mean time-frequency map for one cell for spontaneous and evoked
%upstates
% all signals must be cut 100ms before stimulation and 1 sec after

function [m_TFVmSpont,m_TFEEGSpont,m_TFVmEvok,m_TFEEGEvok,v_TimeAxis]=...
    TFabs1sec_spontevok(st_info,b_spont,b_evok,s_TFreso,b_norm)

    params = getParams();
    s_baseline = params.s_baseline;
    
    b_log = false;
    
    if nargin == 1
        b_spont = true;
        b_evok = true;
    end
    if nargin == 1 || nargin == 3
        s_TFreso = params.s_TFreso;
    end
    if nargin == 1|| nargin == 3 || nargin == 4
        b_norm = true;
    end
    
    s_normBin = 1;
    if b_norm
        s_freqMin = 8;
    else s_freqMin = params.s_freqMin;
    end
    
    %extract upstates
    if st_info(2).header.sampleinterval ~= 360
        [p,q] = rat(st_info(1).header.sampleinterval/60,0.001);
        st_info(1).data = resample(st_info(1).data,p,q);
        st_info(1).header.sampleinterval = 60;
        st_info(1).header.Sampling = 1.666666666666667e+04;
        [p,q] = rat(st_info(2).header.sampleinterval/360,0.001);
        st_info(2).data = resample(st_info(2).data,p,q);
        st_info(2).header.sampleinterval = 360;
        st_info(2).header.Sampling = 2.777777777777778e+03;
    end
    
    [m_spont,m_evok,st_info] = spont_evok_extract(st_info,b_spont,b_evok,true);
    s_sampleInterval = st_info(1).header.sampleinterval
    s_samplingfreq = st_info(1).header.Sampling

      
    m_TFVmSpont = [];
    m_TFEEGSpont = [];
    s_count = 0;
    if b_spont
        
        for s_which = 1:length(m_spont)
            v_unitVmSpont = m_spont{s_which,1};
            v_unitEEGSpont = m_spont{s_which,2};
            [b_valid,v_unitVmSpont] = validity_test(v_unitVmSpont,...
                s_sampleInterval);
            if b_valid
        %% Compute time-frequency map
                s_count = s_count + 1;
                
                %time-frequency map
                s_cut = floor(s_baseline/st_info(1).header.sampleinterval*1000);
                [m_vmSpont,v_TimeAxis,v_FreqAxis] = compute_TF(v_unitVmSpont,...
                    st_info(1).header.Sampling,s_freqMin,100,s_TFreso);
                m_EEGSpontL = compute_TF(v_unitEEGSpont(1:s_cut),...
                    st_info(2).header.Sampling,s_freqMin,100,s_TFreso);
                m_EEGSpontG = compute_TF(v_unitEEGSpont(s_cut+1:end),...
                    st_info(2).header.Sampling,s_freqMin,100,s_TFreso);
                
                m_EEGSpont = [m_EEGSpontL,m_EEGSpontG];
                
                if b_log
                    m_vmSpont = 20*log10(m_vmSpont)
                    m_EEGSpont = 20*log10(m_EEGSpont);
                end
                
                v_TimeAxis = v_TimeAxis(1:length(m_vmSpont));
                
                if b_norm
                    %Normalize amplitude
                    for index = 0:floor(length(v_FreqAxis)/s_normBin)-1
                        m_vmSpont(s_normBin*index+1:s_normBin*(index+1),:) = ...
                            m_vmSpont(s_normBin*index+1:s_normBin*(index+1),:)/...
                            sum(sum(m_vmSpont(s_normBin*index+1:s_normBin*(index+1),:)));
                       m_EEGSpont(s_normBin*index+1:s_normBin*(index+1),:) = ...
                            m_EEGSpont(s_normBin*index+1:s_normBin*(index+1),:)/...
                            sum(sum(m_EEGSpont(s_normBin*index+1:s_normBin*(index+1),:)));
                    end
                end
                
        %% compute average map %%
                if isempty(m_TFVmSpont)
                    m_TFVmSpont = m_vmSpont;
                    m_TFEEGSpont = m_EEGSpont;
                else m_TFVmSpont = m_TFVmSpont + m_vmSpont;
                     m_TFEEGSpont = m_TFEEGSpont + m_EEGSpont;
                end
                    
                
%% Visualize                
%                  figure(14)
%                  subplot(2,1,1)
%                  v_clim = [prctile(m_vmSpont(:),1),prctile(m_vmSpont(:),100)];
%                  f_ImageArray(m_vmSpont,v_TimeAxis,v_FreqAxis,'colormap','jet','limits',v_clim);
%                  xlabel('Temps (s)'); ylabel('Fréquence (Hz)'); title('Intra Spontané')
%                  subplot(2,1,2)
%                  plot(v_TimeAxis,v_unitVmSpont)
%                  saveas(gcf,['/home/eloise/Stages/Rythm-ICM/Matlab/1176-trials/tfIntraSpont_',int2str(s_which),'.png']);

%                  figure(15)
%                  subplot(2,1,1)
%                  v_clim = [prctile(m_EEGSpont(:),1),prctile(m_EEGSpont(:),100)];
%                  f_ImageArray(m_EEGSpont,v_TimeAxis,v_FreqAxis,'colormap','jet','limits',v_clim);
%                  xlabel('Temps (s)'); ylabel('Fréquence (Hz)'); title('EEG Spontané')
%                  subplot(2,1,2)
%                  plot(v_TimeAxis,v_unitEEGSpont)
%                  if b_norm
%                      saveas(gcf,['/home/eloise/Stages/Rythm-ICM/Matlab/1176-trials/tfEEGSpontNorm_',int2str(s_which),'.png']);
%                  else
%                      saveas(gcf,['/home/eloise/Stages/Rythm-ICM/Matlab/1176-trials/tfEEGSpont_',int2str(s_which),'.png']);
%                  end
                                  
            end
            
            
        end
      
        
        m_TFVmSpont = m_TFVmSpont/s_count;
        m_TFEEGSpont = m_TFEEGSpont/s_count;
        %% normalize average map
%         for index = 1:length(v_FreqAxis)
%             m_TFVmSpont(index,:) = m_TFVmSpont(index,:)/sum(m_TFVmSpont(index,:));
%             m_TFEEGSpont(index,:) = m_TFEEGSpont(index,:)/sum(m_TFEEGSpont(index,:));
%         end
    end
%% Évoqué    
    m_TFVmEvok = [];
    m_TFEEGEvok = [];
    s_count = 0;
    if b_evok
        for s_which = 1:length(m_evok)
            v_unitVmEvok = m_evok{s_which,1}; 
            v_unitEEGEvok = m_evok{s_which,2};
            [b_valid,v_unitVmEvok] = validity_test(v_unitVmEvok,s_sampleInterval);
            if b_valid
         %% Compute time_frequency map
                s_count = s_count + 1;
                %cut steep slope if necessary
                [m_vmEvok,v_TimeAxis] = TF_evok_cut_all(v_unitVmEvok,'Vm',...
                    s_sampleInterval,s_samplingfreq,s_TFreso,s_freqMin);
                m_EEGEvok = TF_evok_cut_all(v_unitEEGEvok,'EEG',s_sampleInterval,...
                    s_samplingfreq, s_TFreso,s_freqMin);
                v_FreqAxis = linspace(100,s_freqMin,s_TFreso);
                if b_log
                    m_vmEvok = 20*log10(m_vmEvok);
                    m_EEGEvok = 20*log10(m_EEGEvok);
                end        
                if b_norm
                    %Normalize
                   for index = 0:floor(length(v_FreqAxis)/s_normBin)-1
                            m_vmEvok(s_normBin*index+1:s_normBin*(index+1),:) = ...
                                m_vmEvok(s_normBin*index+1:s_normBin*(index+1),:)/...
                                sum(sum(m_vmEvok(s_normBin*index+1:s_normBin*(index+1),:)));
                           m_EEGEvok(s_normBin*index+1:s_normBin*(index+1),:) = ...
                                m_EEGEvok(s_normBin*index+1:s_normBin*(index+1),:)/...
                                sum(sum(m_EEGEvok(s_normBin*index+1:s_normBin*(index+1),:)));
                    end
                end   
    
%% Visualize                
%              figure(14)
%              subplot(2,1,1)
%              v_clim = [prctile(m_vmEvok(:),1),prctile(m_vmEvok(:),100)];
%              f_ImageArray(m_vmEvok,v_TimeAxis,v_FreqAxis,'colormap','jet','limits',v_clim);
%              xlabel('Temps (s)'); ylabel('Fréquence (Hz)'); title('Intra Évoqué');
%              subplot(2,1,2)
%              plot(v_TimeAxis,v_unitVmEvok)
%              saveas(gcf,['/home/eloise/Stages/Rythm-ICM/Matlab/1176-trials/tfIntraEvok_',int2str(s_which),'.png']);
% 
%              figure(15)
%              subplot(2,1,1)
%              v_clim = [prctile(m_EEGEvok(:),1),prctile(m_EEGEvok(:),100)];
%              f_ImageArray(m_EEGEvok,v_TimeAxis,v_FreqAxis,'colormap','jet','limits',v_clim);
%              xlabel('Temps (s)'); ylabel('Fréquence (Hz)'); title('EEG Évoqué');
%              subplot(2,1,2)
%              plot(v_TimeAxis,v_unitEEGEvok)
%              if b_norm
%                  saveas(gcf,['/home/eloise/Stages/Rythm-ICM/Matlab/1176-trials/tfEEGEvokNorm_',int2str(s_which),'.png']);
%              else 
%                  saveas(gcf,['/home/eloise/Stages/Rythm-ICM/Matlab/1176-trials/tfEEGEvok_',int2str(s_which),'.png']);
%              end
             
        %% Average
                if isempty(m_TFVmEvok)
                    m_TFVmEvok = m_vmEvok;
                    m_TFEEGEvok = m_EEGEvok;
                else m_TFVmEvok = m_TFVmEvok + m_vmEvok;
                     m_TFEEGEvok = m_TFEEGEvok + m_EEGEvok;
                end

            end
        end
        
        m_TFVmEvok = m_TFVmEvok/s_count;
        m_TFEEGEvok = m_TFEEGEvok/s_count;
        %% Normalize average map
    %      for index = 1:length(v_FreqAxis)
    %             m_TFVmEvok(index,:) = m_TFVmEvok(index,:)/sum(m_TFVmEvok(index,:));
    %             m_TFEEGEvok(index,:) = m_TFEEGEvok(index,:)/sum(m_TFEEGEvok(index,:));
    %         end
        
    end
end
