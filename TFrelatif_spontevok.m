function [m_TFVmSpont,m_TFEEGSpont,m_TFVmEvok,m_TFEEGEvok]=...
    TFrelatif_spontevok(info,b_spont,b_evok,s_nbBins,s_TFreso,b_norm)

    params = getParams();

    if nargin == 1
        b_spont = true;
        b_evok = true;
    end
    
    if nargin == 1 || nargin == 3
        s_nbBins = params.s_nbBins;
    end
    
    if nargin == 1 || nargin == 3 || nargin == 4
        s_TFreso = params.s_TFreso;
    end
    
    if nargin == 1|| nargin == 3 || nargin == 4 || nargin == 5
        b_norm = true;
    end
    
    if nargin == 2
         error('[TFrelatif_spontevok] - ERROR: Bad number of parameters')   
    end
    
      if b_norm
        s_freqMin = 0;
    else s_freqMin = params.s_freqMin;
      end
    
      
    [m_spont,m_evok,info] = spont_evok_extract(info,b_spont,b_evok);
    s_sampleInterval = info(1).header.sampleinterval;
    s_samplingfreq = info(1).header.Sampling;

    m_TFVmSpont = [];
    m_TFEEGSpont = [];
    m_TFVmEvok = [];
    m_TFEEGEvok = [];
    s_count=0;
    
    
    if b_spont    
        for s_which=1:length(m_spont)
            v_unitVmSpont =m_spont{s_which,1}; 
            v_unitEEGSpont = m_spont{s_which,2};
            [b_valid,v_unitVmSpont] = validity_test(v_unitVmSpont,s_sampleInterval);
            if b_valid
%% Compute time-frequency map                
                s_count = s_count + 1;                
                m_vmSpont = compute_TF(v_unitVmSpont,info(1).header.Sampling,s_freqMin,100,s_TFreso);
                [m_EEGSpont, v_TimeAxis, v_FreqAxis] = ...
                compute_TF(v_unitEEGSpont,info(2).header.Sampling,s_freqMin,100,s_TFreso);

                %log scale
                %m_vmSpont = 20*log10(m_vmSpont(:,floor(end/3):floor(2*end/3)-1));
                %m_EEGSpont = 20*log10(m_EEGSpont(:,floor(end/3):floor(2*end/3)-1));
                
                if b_norm
                    %Normalize
                    for index = 1:length(v_FreqAxis)
                        m_vmSpont(index,:) = m_vmSpont(index,:)/sum(m_vmSpont(index,:));
                        m_EEGSpont(index,:) = m_EEGSpont(index,:)/sum(m_EEGSpont(index,:));
                    end
                end
%% Compute average
                m_redVmSpont = zeros(size(m_vmSpont,1),s_nbBins);
                m_redEEGSpont = zeros(size(m_EEGSpont,1),s_nbBins);
                s_W = size(m_vmSpont,2);
                s_WBin =floor(s_W/s_nbBins);
                v_aux = 1:s_WBin:s_nbBins*s_WBin;
                v_TimeAxis = v_TimeAxis(v_aux);
                for i=1:s_nbBins
                    m_redEEGSpont(:,i) = mean(m_EEGSpont(:,(i-1)*s_WBin +1: i*s_WBin),2);
                    m_redVmSpont(:,i) = mean(m_vmSpont(:,(i-1)*s_WBin +1: i*s_WBin),2);
                end
                if isempty(m_TFVmSpont)
                    m_TFVmSpont =m_redVmSpont;
                else m_TFVmSpont = m_TFVmSpont + m_redVmSpont;
                end
                 if isempty(m_TFEEGSpont)
                    m_TFEEGSpont =m_redEEGSpont;
                else m_TFEEGSpont = m_TFEEGSpont + m_redEEGSpont;
                 end
%% Visualize
%                 time = 0:s_sampleInterval/1000000:s_sampleInterval/1000000*(length(v_unitVmSpont)-1);
%                 figure(12)
%                 subplot(1,3,1);
%                     f_ImageArray(m_redVmSpont,v_TimeAxis, v_FreqAxis,'colormap','jet');
%                     title('Intra');xlabel('temps (s)');ylabel('Fr�quence (Hz)');
%                 subplot(1,3,2);
%                     f_ImageArray(m_redEEGSpont,v_TimeAxis, v_FreqAxis,'colormap','jet');
%                     title('EEG');xlabel('temps (s)');ylabel('Fr�quence (Hz)');
%                 subplot(2,3,3);
%                     plot(time,v_unitVmSpont); title('Intra');xlabel('temps (s');ylabel('Potentiel (mV)');
%                 subplot(2,3,6);
%                     plot(time,v_unitEEGSpont); title('EEG');xlabel('temps (s');ylabel('Potentiel (muV)');
%     
%                 figure(13)
%                 subplot(2,2,1);
%                     f_ImageArray(m_TFVmSpont,v_TimeAxis, v_FreqAxis,'colormap','jet');
%                     title('Vm');xlabel('temps (s)');ylabel('Fr�quence (Hz)');
%                 subplot(2,2,2);
%                     f_ImageArray(m_TFEEGSpont,v_TimeAxis, v_FreqAxis,'colormap','jet');
%                     title('EEG');xlabel('temps (s)');ylabel('Fr�quence (Hz)');
%                 subplot(3,2,5)
%                     plot(time,v_unitVmSpont);
%                     hold on;
%                     title('Intra');xlabel('temps (s)');ylabel('Potentiel (mV)');
%                 subplot(3,2,6)
%                     plot(time,v_unitEEGSpont);
%                     title('EEG');xlabel('temps (S)');ylabel('Potentiel (mV)');
%                     hold on;
                
            end
        end
        %% Normalize average map
%         for index = 1:length(v_FreqAxis)
%             m_TFVmSpont(index,:) = m_TFVmSpont(index,:)/sum(m_TFVmSpont(index,:));
%             m_TFEEGSpont(index,:) = m_TFEEGSpont(index,:)/sum(m_TFVmSpont(index,:));
%         end
    end
 
    if b_evok
        s_count=0;
        for s_which=1:length(m_evok)
            v_unitVmEvok = m_evok{s_which,1}; 
            v_unitEEGEvok = m_evok{s_which,2};
            %b_valid = true;
            [b_valid,v_unitVmEvok] = validity_test(v_unitVmEvok,s_sampleInterval);
            if b_valid
%% Compute time-frequency map

                s_count=s_count+1;
                m_vmEvok = compute_TF(v_unitVmEvok,info(1).header.Sampling,s_freqMin,100,s_TFreso);
                [m_EEGEvok, v_TimeAxis, v_FreqAxis] = ...
                fcompute_TF(v_unitEEGEvok,info(2).header.Sampling,s_freqMin,100,s_TFreso);
                %log scale
                %m_vmEvok = 20*log10(m_vmEvok);
                %m_EEGEvok = 20*log10(m_EEGEvok);
                
                if b_norm
                    %Normalize
                    for index = 1:length(v_FreqAxis)
                        m_vmEvok(index,:) = m_vmEvok(index,:)/sum(m_vmEvok(index,:));
                        m_EEGEvok(index,:) = m_EEGEvok(index,:)/sum(m_EEGEvok(index,:));
                    end
                end
%% Compute average
                m_redVmevok = zeros(size(m_vmEvok,1),s_nbBins);
                m_redEEGevok = zeros(size(m_EEGEvok,1),s_nbBins);
                s_W = size(m_vmEvok,2);
                s_WBin =floor(s_W/s_nbBins);
                v_aux = 1:s_WBin:s_nbBins*s_WBin;
                v_TimeAxis = v_TimeAxis(v_aux);
                for i=1:s_nbBins
                    m_redEEGevok(:,i) = mean(m_EEGEvok(:,(i-1)*s_WBin +1: i*s_WBin),2);
                    m_redVmevok(:,i) = mean(m_vmEvok(:,(i-1)*s_WBin +1: i*s_WBin),2);
                end
                if isempty(m_TFVmEvok)
                    m_TFVmEvok =m_redVmevok;
                else m_TFVmEvok = m_TFVmEvok + m_redVmevok;
                end
                 if isempty(m_TFEEGEvok)
                    m_TFEEGEvok =m_redEEGevok;
                else m_TFEEGEvok = m_TFEEGEvok + m_redEEGevok;
                 end
%% Visualize
                time = 0:s_sampleInterval/1000000:s_sampleInterval/1000000*(length(v_unitVmEvok)-1); 
%                 figure(14)
%                 subplot(1,3,1);
%                     f_ImageArray(m_redVmevok,v_TimeAxis, v_FreqAxis,'colormap','jet');
%                     title('Intra');xlabel('temps (s)');ylabel('Fr�quence (Hz)');
%                 subplot(1,3,2);
%                     f_ImageArray(m_redEEGevok,v_TimeAxis, v_FreqAxis,'colormap','jet');
%                     title('EEG');xlabel('temps (s)');ylabel('Fr�quence (Hz)');
%                 subplot(2,3,3);
%                     plot(time,v_unitVmEvok); title('Intra');xlabel('temps (s');ylabel('Potentiel (mV)');
%                 subplot(2,3,6);
%                     plot(time,v_unitEEGEvok); title('EEG');xlabel('temps (s');ylabel('Potentiel (muV)');
% 
%                 figure(15)
%                 subplot(2,2,1);
%                     f_ImageArray(m_TFVmEvok,v_TimeAxis, v_FreqAxis,'colormap','jet');
%                     title('Vm');xlabel('temps (s)');ylabel('Fr�quence (Hz)');
%                 subplot(2,2,2);
%                     f_ImageArray(m_TFEEGEvok,v_TimeAxis, v_FreqAxis,'colormap','jet');
%                     title('EEG');xlabel('temps (s)');ylabel('Fr�quence (Hz)');
%                 subplot(3,2,5)
%                     plot(time,v_unitVmEvok);
%                     hold on;
%                     title('Intra');xlabel('temps (s)');ylabel('Potentiel (mV)');
%                 subplot(3,2,6)
%                     plot(time,v_unitEEGEvok);
%                     title('EEG');xlabel('temps (s)');ylabel('Potentiel (mV)');
%                     hold on;
% 
%            keyboard;

            end
        end
        %% Normalize average map
%         for index = 1:length(v_FreqAxis)
%             m_TFVmEvok(index,:) = m_TFVmEvok(index,:)/sum(m_TFVmEvok(index,:));
%             m_TFEEGEvok(index,:) = m_TFEEGEvok(index,:)/sum(m_TFEEGEvok(index,:));
%         end
    end
end
