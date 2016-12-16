%Compute time frequency map with given frequency band and resolution


function [m_TF,v_TimeAxis,v_FreqAxis] = compute_TF(v_signal,s_sampling,...
    s_freqmin,s_freqmax,s_TFreso)

        %mirror before wavelet transform to avoid edge effect
        v_auxG = 2*v_signal(1)-v_signal(end:-1:1);
        v_auxR = 2*v_signal(end)-v_signal(end:-1:1);
        v_signal3 = [v_auxG,v_signal,v_auxR];
        
        [m_TF, v_TimeAxis, v_FreqAxis] = ...
                    f_GaborTransformWait(v_signal3,...
                    s_sampling,s_freqmin,s_freqmax,s_TFreso);
                
         m_TF = m_TF(:,floor(end/3):floor(2*end/3)-1);
         v_TimeAxis = v_TimeAxis(1:length(m_TF));


end
