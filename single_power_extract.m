% Extract power in given frequency band on given time frequency map

function [v_Spont,v_Evok] = ...
    single_power_extract(m_TFAbsVmSpont,m_TFAbsVmEvok,v_FreqAxis,v_freqBand)

    if numel(v_freqBand)<2
         error('[single_power_extract] - ERROR: invalid frequency band')
    end

    s_freqMin = v_freqBand(1);
    s_freqMax = v_freqBand(2);
    v_FreqAxis = v_FreqAxis(end:-1:1); 
    m_TFAbsVmSpont = m_TFAbsVmSpont(end:-1:1,:);
    m_TFAbsVmEvok = m_TFAbsVmEvok(end:-1:1,:);
    index = 1;
    v_Spont = zeros(1,size(m_TFAbsVmSpont,2));
    v_Evok = zeros(1,size(m_TFAbsVmEvok,2));

    while v_FreqAxis(index)<s_freqMin
       index = index + 1;
    end

    start = index;
    while v_FreqAxis(index)<s_freqMax
        v_Spont = v_Spont + m_TFAbsVmSpont(index,:);
        v_Evok = v_Evok + m_TFAbsVmEvok(index,:);
        index = index + 1;
    end
    v_Spont = v_Spont/(index - start);
    v_Evok = v_Evok/(index - start);
end

