%Filters one upstate and the corresponding EEG on frequency interval :
%spindle, gamma Low(30-60Hz), gamma Middle(60-80Hz) and gamma High(80-100Hz)

function [v_FiltSpVm,v_FiltGamLoVm,v_FiltGamMidVm, v_FiltGamHiVm,...
    v_FiltSpEEG,v_FiltGamLoEEG,v_FiltGamMidEEG,v_FiltGamHiEEG] =...
    filtre_upstate(v_Vm,v_EEG,s_sampling_freq)


    %mirror before filtering to avoid edge effet
    v_aux_vm_G = 2*v_Vm(1)-v_Vm(end:-1:1);
    v_aux_vm_R = 2*v_Vm(end)-v_Vm(end:-1:1);
    v_filt_vm = [v_aux_vm_G,v_Vm,v_aux_vm_R];
    v_aux_EEG_G = 2*v_EEG(1)-v_EEG(end:-1:1);
    v_aux_EEG_R = 2*v_EEG(end)-v_EEG(end:-1:1);
    v_filt_EEG = [v_aux_EEG_G,v_EEG,v_aux_EEG_R];

    %Filtering
    v_FiltSpindle	= f_DesignIIRfilter(s_sampling_freq,[9 16],[8 17],[1 100]); 
    v_FiltGammaLo	= f_DesignIIRfilter(s_sampling_freq,[30 60],[29 61],[1,50]);
    v_FiltGammaMid	= f_DesignIIRfilter(s_sampling_freq,[60 80],[59 81],[1,50]);
    v_FiltGammaHi	= f_DesignIIRfilter(s_sampling_freq,[80 100],[79 101],[1,50]);
    % Intra
    v_FiltSpVm	= f_FilterIIR(v_filt_vm,v_FiltSpindle);
    v_FiltSpVm = v_FiltSpVm(floor(end/3):floor(2*end/3));
    v_FiltGamLoVm = f_FilterIIR(v_filt_vm,v_FiltGammaLo);
    v_FiltGamLoVm = v_FiltGamLoVm(floor(end/3):floor(2*end/3));
    v_FiltGamMidVm = f_FilterIIR(v_filt_vm,v_FiltGammaMid); 
    v_FiltGamMidVm = v_FiltGamMidVm(floor(end/3):floor(2*end/3));
    v_FiltGamHiVm = f_FilterIIR(v_filt_vm,v_FiltGammaHi); 
    v_FiltGamHiVm = v_FiltGamHiVm(floor(end/3):floor(2*end/3));
    %EEG 
    v_FiltSpEEG	= f_FilterIIR(v_filt_EEG,v_FiltSpindle);
    v_FiltSpEEG = v_FiltSpEEG(floor(end/3):floor(2*end/3));
    v_FiltGamLoEEG = f_FilterIIR(v_filt_EEG,v_FiltGammaLo);
    v_FiltGamLoEEG = v_FiltGamLoEEG(floor(end/3):floor(2*end/3));
    v_FiltGamMidEEG = f_FilterIIR(v_filt_EEG,v_FiltGammaMid);
    v_FiltGamMidEEG = v_FiltGamMidEEG(floor(end/3):floor(2*end/3));
    v_FiltGamHiEEG = f_FilterIIR(v_filt_EEG,v_FiltGammaHi);
    v_FiltGamHiEEG = v_FiltGamHiEEG(floor(end/3):floor(2*end/3));

end