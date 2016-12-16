% Cut evoked signal where the depolarization slope is so high that it will
% cause an artefact on the time frequency map, and then compute left and
% right frequency maps

function [m_totalEvok,v_totalTimeAxis] = TF_evok_cut(...
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
    
    %% Compute TF map without cutting slope

    %time-frequency map 
    [m_Evok, v_totalTimeAxis, v_FreqAxis] = compute_TF(v_Evok...
                    ,s_samplingfreq,s_freqMin,100,s_TFreso);

%% Define depolarization slope interval
  
    %find abrupt slope of more than one sample interval
    %we know that it should be between 100 ms and 200 ms
    s_limDown = floor(s_baseline/s_sampleinterval*1000);
    s_limUp = floor((s_baseline+100)/s_sampleinterval*1000);
    switch signaltype
        case 'Vm'
            v_diff = (v_evok(2:end) - v_evok(1:end-1))/360*1000;
            v_as = find(abs(v_diff(s_limDown:s_limUp))>1.5)+s_limDown;
        case 'EEG'
            v_evoknorm = (v_evok - mean(v_evok))/std(v_evok);
            v_diff = (v_evoknorm(2:end) - v_evoknorm(1:end-1))/360*1000;
            v_as = find(abs(v_diff(s_limDown:s_limUp))>0.5)+s_limDown;
    end
    
    % if isolated sample remove
    v_asnew = v_as;
    for idx = 1:length(v_as)
        if (idx==1 || v_as(idx-1)<v_as(idx)-1) && (idx==length(v_as) || v_as(idx+1)>v_as(idx)+1)
            v_asnew(idx)=0;
        end
    end
    v_as = v_as(find(v_asnew));
    
    %dont confuse with following slope
    v_chge = find(v_as(2:end)-v_as(1:end-1)>2); 
    if ~isempty(v_chge)
        v_as = v_as(1:v_chge(1));
    end
    
    if length(v_as) >3 % otherwise no artefact
        v_as(1) = v_as(1) - 15; %security interval
        v_as(end) = v_as(end) + 15;
        v_leftEvok = v_evok(1:v_as(1));
        v_rightEvok = v_evok(v_as(end)+1:end);
        %visualize cut
%         v_test = [v_leftEvok,zeros(1,v_as(end)-v_as(1)),v_rightEvok];
%         figure(3)
%         plot(v_evok)
%         hold on;
%         plot(v_test)
%         hold off;

%% Compute left and right TF maps
        %mirror before wavelet transform to avoid edge effect
        v_auxG = 2*v_leftEvok(1)-v_leftEvok(end:-1:1);
        v_auxR = 2*v_leftEvok(end)-v_leftEvok(end:-1:1);
        v_leftEvok3 = [v_auxG,v_leftEvok,v_auxR];
        %time-frequency map
        [m_leftEvok,v_TimeAxisleft,v_FreqAxis]=...
            f_GaborTransformWait(v_leftEvok3,s_samplingfreq,s_freqMin,100,s_TFreso);
        %cut mirror
        m_leftEvok = m_leftEvok(:,floor(end/3):floor(2*end/3)-1);

        %mirror before wavelet transform to avoid edge effect
        v_auxG = 2*v_rightEvok(1)-v_rightEvok(end:-1:1);
        v_auxR = 2*v_rightEvok(end)-v_rightEvok(end:-1:1);
        v_rightEvok3 = [v_auxG,v_rightEvok,v_auxR];
        %time-frequency map
        [m_rightEvok,v_TimeAxisright,~]=...
            f_GaborTransformWait(v_rightEvok3,s_samplingfreq,s_freqMin,100,s_TFreso);
        %cut mirror
        m_rightEvok = m_rightEvok(:,floor(end/3):floor(2*end/3)-1);
        %coeff = max(max(m_rightEvok))/max(max(m_leftEvok));
        m_totalEvok = [m_leftEvok,ones(length(v_FreqAxis),v_as(end)-v_as(1)),m_rightEvok];
    else
        m_totalEvok = m_Evok;
    end
    
%% Visualize difference
%     figure(4)
%     subplot(1,2,1)
%         f_ImageArray(m_Evok,v_totalTimeAxis,v_FreqAxis,'colormap','jet');
%     subplot(1,2,2)
%         f_ImageArray(m_totalEvok,v_totalTimeAxis,v_FreqAxis,'colormap','jet');
%         keyboard;
        
end