%Computes power spectra in gamma and spindle frequencies from a given
%upstate

function [v_powSpVm,v_powGamVm,v_powSpEEG,v_powGamEEG] = ...
    filtre_spectre(v_Vm,v_EEG,s_sampling_freq)
    
%% Filtering
    %mirror on both sides before filtering to avoid edge effect
    v_auxVmG = 2*v_Vm(1)-v_Vm(end:-1:1);
    v_auxVmR = 2*v_Vm(end)-v_Vm(end:-1:1);
    v_unitVm = [v_auxVmG,v_Vm,v_auxVmR];
    v_auxEEGG = 2*v_EEG(1)-v_EEG(end:-1:1);
    v_auxEEGR = 2*v_EEG(end)-v_EEG(end:-1:1);
    v_unitEEG = [v_auxEEGG,v_EEG,v_auxEEGR];

    v_FiltSpindle	= f_DesignIIRfilter(s_sampling_freq,[9 16],[8 17],[1 100]); 
    v_FiltGamma	= f_DesignIIRfilter(s_sampling_freq,[30 100],[29 101],[1,100]);
    % Intra
    v_FiltSpVm	= f_FilterIIR(v_unitVm,v_FiltSpindle);
    v_FiltSpVm = v_FiltSpVm(floor(end/3):floor(2*end/3));
    v_FiltGamVm = f_FilterIIR(v_unitVm,v_FiltGamma); 
    v_FiltGamVm = v_FiltGamVm(floor(end/3):floor(2*end/3));
    %EEG 
    v_FiltSpEEG	= f_FilterIIR(v_unitEEG,v_FiltSpindle);
    v_FiltSpEEG = v_FiltSpEEG(floor(end/3):floor(2*end/3));
    v_FiltGamEEG = f_FilterIIR(v_unitEEG,v_FiltGamma);
    v_FiltGamEEG = v_FiltGamEEG(floor(end/3):floor(2*end/3));

    time=0:0.36:(length(v_FiltSpEEG)-1)*0.36;
    figure(11)
    subplot(2,2,1)
        plot(time,v_FiltSpVm);
        title('Spindle intra')
        xlabel('Temps (ms)')
        ylabel('Amplitude (mV)')
    subplot(2,2,2)
        plot(time,v_FiltSpEEG)
        title('Spindle EEG')
        xlabel('Temps (ms)')
        ylabel('Amplitude (muV)')
    subplot(2,2,3)
        plot(time,v_FiltGamVm)
        xlabel('Temps (ms)')
        ylabel('Amplitude (mV)')
        title('Gamma intra')
    subplot(2,2,4)
        plot(time,v_FiltGamEEG)
        xlabel('Temps (ms)')
        ylabel('Amplitude (muV)')
        title('Gamma EEG')

%% Power spectrum %%
    NFFT = 2^nextpow2(300/0.36);
    f = s_sampling_freq/2*linspace(0,1,NFFT/2+1);

    %Vm
    v_fft_spindle_vm = fft(v_FiltSpVm,NFFT);
    v_fft_spindle_vm = v_fft_spindle_vm(1:floor(end/2)+1);
    v_fft_gamma_vm = fft(v_FiltGamVm,NFFT);
    v_fft_gamma_vm = v_fft_gamma_vm(1:floor(end/2)+1);
    v_powSpVm = v_fft_spindle_vm.*conj(v_fft_spindle_vm);
    v_powGamVm = v_fft_gamma_vm.*conj(v_fft_gamma_vm);


    %EEG
    v_fft_spindle_EEG = fft(v_FiltSpEEG,NFFT);
    v_fft_spindle_EEG = v_fft_spindle_EEG(1:floor(end/2)+1);
    v_fft_gamma_EEG = fft(v_FiltGamEEG,NFFT);
    v_fft_gamma_EEG = v_fft_gamma_EEG(1:floor(end/2)+1);
    v_powSpEEG = v_fft_spindle_EEG.*conj(v_fft_spindle_EEG);
    v_powGamEEG = v_fft_gamma_EEG.*conj(v_fft_gamma_EEG);

    figure(3)
    subplot(2,2,1);
    plot(f,v_powSpVm);
    title('Spectre puissance Spindle vm');
    subplot(2,2,2);
    plot(f,v_powGamVm);
    title('Spectre puissance Gamma vm');
    subplot(2,2,3);
    plot(f,v_powSpEEG);
    title('Spectre puissance spindle EEG');
    subplot(2,2,4);
    plot(f,v_powGamEEG);
    title('Spectre puissance gamma EEG');
    

end