%open info file to extract Vm and EEG and subsample to the same Sampling

function [st_info, v_Vm,v_EEG]=open_subsamp_VmEEG(st_info)    
    
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
    
    v_VmOrigin = st_info(1).data;
    v_EEG = st_info(2).data;
    
    %setting intra and EEG to the same sampling rate
    s_subsampRate = st_info(2).header.sampleinterval/st_info(1).header.sampleinterval;
    v_index = 1:s_subsampRate:length(v_VmOrigin);
    st_info(1).header.Sampling = st_info(1).header.Sampling/s_subsampRate;
    st_info(1).header.sampleinterval = st_info(1).header.sampleinterval*s_subsampRate;
    v_Vm = v_VmOrigin(v_index);
    %assert equal length
    s_min = min(length(v_EEG),length(v_Vm));
    v_Vm = v_Vm(1:s_min); 
    v_EEG = v_EEG(1:s_min);
    st_info(1).data = v_Vm;
    st_info(2).data = v_EEG;
    assert(length(v_Vm)==length(v_EEG));
    
    
    
    
end
