%correlations
info = openSpike2('period1.smr');
%découpe les upstates dans EEG et intra
[m_upstates,info] = upstate_analysis(info);


Si = info(1).header.sampleinterval;
Fe = info(1).header.Sampling;
nb_windows = 50;


for s_which=1:length(m_upstates)
%sélection d'un up de travail
v_vm_upstates= m_upstates{s_which,1};
v_EEG_upstates = m_upstates{s_which,2};
b_len = length(v_vm_upstates)*Si/1000 > 300;
v_diff = (v_vm_upstates(2:end) - v_vm_upstates(1:end-1))/360*1000;
b_pa = sum(v_diff >10) <0.5;


if b_len && b_pa
  
    [v_FiltSpindle_vm,v_FiltGammaLo_vm,v_FiltGammaMid_vm, v_FiltGammaHi_vm,...
    v_FiltSpindle_EEG,v_FiltGammaLo_EEG,v_FiltGammaMid_EEG,v_FiltGammaHi_EEG] =...
    filtre_upstate(v_vm_upstates,v_EEG_upstates,info(1).header.Sampling);
   size_wd = 2*floor(length(v_FiltSpindle_vm)/nb_windows);

     x_sp = v_FiltSpindle_vm - mean(v_FiltSpindle_vm);
     y_sp = v_FiltSpindle_EEG - mean(v_FiltSpindle_EEG);
     x_gamLo = v_FiltGammaLo_vm - mean(v_FiltGammaLo_vm);
     y_gamLo = v_FiltGammaLo_EEG - mean(v_FiltGammaLo_EEG);
     x_gamMid = v_FiltGammaMid_vm - mean(v_FiltGammaMid_vm);
     y_gamMid = v_FiltGammaMid_EEG - mean(v_FiltGammaMid_EEG);
     x_gamHi = v_FiltGammaHi_vm - mean(v_FiltGammaHi_vm);
     y_gamHi = v_FiltGammaHi_EEG - mean(v_FiltGammaHi_EEG);
     v_xcorr_spindle=zeros(nb_windows-1);
     v_xcorr_gamLo=zeros(nb_windows-1);
     v_xcorr_gamMid=zeros(nb_windows-1);
     v_xcorr_gamHi=zeros(nb_windows-1);
     
     for i =1:(nb_windows-1)
        x_sp_tr = x_sp((i-1)*size_wd/2+1:(i-1)*size_wd/2+size_wd);
        y_sp_tr = y_sp((i-1)*size_wd/2+1:(i-1)*size_wd/2+size_wd);
        x_gamLo_tr = x_gamLo((i-1)*size_wd/2+1:(i-1)*size_wd/2+size_wd);
        y_gamLo_tr = y_gamLo((i-1)*size_wd/2+1:(i-1)*size_wd/2+size_wd);
        x_gamMid_tr = x_gamMid((i-1)*size_wd/2+1:(i-1)*size_wd/2+size_wd);
        y_gamMid_tr = y_gamMid((i-1)*size_wd/2+1:(i-1)*size_wd/2+size_wd);
        x_gamHi_tr = x_gamHi((i-1)*size_wd/2+1:(i-1)*size_wd/2+size_wd);
        y_gamHi_tr = y_gamHi((i-1)*size_wd/2+1:(i-1)*size_wd/2+size_wd);
        
        [v_xcorr_spindle_tr,lag_sp] = xcorr(x_sp_tr,y_sp_tr, 'coeff');
        [v_xcorr_gamLo_tr,lag_gamLo] = xcorr(x_gamLo_tr,y_gamLo_tr, 'coeff');
        [v_xcorr_gamMid_tr,lag_gamMid] = xcorr(x_gamMid_tr,y_gamMid_tr, 'coeff');
        [v_xcorr_gamHi_tr,lag_gamHi] = xcorr(x_gamHi_tr,y_gamHi_tr, 'coeff');
        v_xcorr_spindle(i)=max(v_xcorr_spindle_tr);
        v_xcorr_gamLo(i)=max(v_xcorr_gamLo_tr);
        v_xcorr_gamMid(i)=max(v_xcorr_gamMid_tr);
        v_xcorr_gamHi(i)=max(v_xcorr_gamHi_tr);
       
     end   
    
figure(51)
subplot(2,2,1)
plot(v_xcorr_spindle);
title('Spinlde');
subplot(2,2,2)
plot(v_xcorr_gamLo);
title('Gamma 30-60');
subplot(2,2,3)
plot(v_xcorr_gamMid);
title('Gamma 60-80');
subplot(2,2,4)
plot(v_xcorr_gamHi);
title('Gamma 60-80');

figure(52)
subplot(2,1,1)
plot(v_vm_upstates)
subplot(2,1,2)
plot(v_EEG_upstates);
keyboard;
else continue
end

end

