%Cut a given fixed interval on the signal and then compute left and right
%time frequency maps

function [m_totalEvok,v_totalTimeAxis] = TF_evok_cut_all(...
    v_evok,signaltype,s_sampleinterval,s_samplingfreq,s_TFreso,s_freqMin)

    params = getParams();
    
    if nargin <4
        error('[ERROR] TF_evok_cut - invalid number of arguments');
    end
    
    if nargin == 4
        s_TFreso = params.s_TFreso;
    end
    
    if nargin == 4 || nargin == 5
        s_freqMin = 0;
    end
    
    s_baseline = params.s_baseline;
    s_cutWidth = floor(40/s_sampleinterval*1000);
    
    %% Compute TF map without cutting slope
    
    %time-frequency map 
    [m_Evok, v_totalTimeAxis, v_FreqAxis] = compute_TF(v_evok...
                    ,s_samplingfreq,s_freqMin,100,s_TFreso);
  
   %% define left and right signal
    s_stimTime = floor(s_baseline/s_sampleinterval*1000);
    v_leftEvok = v_evok(1:s_stimTime);
    v_rightEvok = v_evok(s_stimTime + s_cutWidth+1:end);
    

%% Compute left and right TF maps
    
        %time-frequency map
        [m_leftEvok,~,v_FreqAxis]=...
            compute_TF(v_leftEvok,s_samplingfreq,s_freqMin,100,s_TFreso);
        %time-frequency map
        m_rightEvok = compute_TF(v_rightEvok,...
            s_samplingfreq,s_freqMin,100,s_TFreso);
        %coeff = max(max(m_rightEvok))/max(max(m_leftEvok));
        m_totalEvok = [m_leftEvok,ones(length(v_FreqAxis),s_cutWidth),m_rightEvok];
        size(m_totalEvok)
        length(v_totalTimeAxis)
   
%% Visualize difference
%     figure(4)
%     subplot(1,2,1)
%         f_ImageArray(m_Evok,v_totalTimeAxis,v_FreqAxis,'colormap','jet');
%     subplot(1,2,2)
%         f_ImageArray(m_totalEvok,v_totalTimeAxis,v_FreqAxis,'colormap','jet');
%         keyboard;
        
end