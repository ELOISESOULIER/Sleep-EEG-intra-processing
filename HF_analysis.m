%Computes medium spectrum on the whole file
s_nbPeriods = 4;

%initialization
[v_meanSpVm,v_meanGamVm,v_meanSpEEG,v_meanGamEEG]= mean_spec('period1.smr');

for s_period=2:s_nbPeriods
    
    file = strcat('period', int2str(s_period),'.smr');
    [v_meanSpVmTr, v_meanGamVmTr, v_meanSpEEGTr,v_meanGamEEGTr]...
        = mean_spec(file);
    v_meanSpVm = v_meanSpVm + v_meanSpVmTr;
    v_meanGamVm = v_meanGamVm + v_meanGamVmTr;
    v_meanSpEEG = v_meanSpEEG + v_meanSpEEGTr;
    v_meanGamEEG = v_meanGamEEG + v_meanGamEEGTr;
    
end

%% Mean, total power and dominant frequency %%
v_meanSpVm = v_meanSpVm/s_nbPeriods;
[~,i] = max(v_meanSpVm);
s_fdSpVm=2*i;
v_spVmCum= cumsum(v_meanSpVm);
s_powSpVm = max(v_spVmCum);

v_meanGamVm = v_meanGamVm/s_nbPeriods;
[~,i]=max(v_meanGamVm);
s_fdGamVm=2*i;
v_gamVmCum= cumsum(v_meanGamVm);
s_powGamVm = max(v_gamVmCum);

v_meanSpEEG = v_meanSpEEG/s_nbPeriods;
[~,i]=max(v_meanSpEEG);
s_fdSpEEG=2*i;
v_spEEGCum= cumsum(v_meanSpEEG);
s_powSpEEG = max(v_spEEGCum);

v_meanGamEEG = v_meanGamEEG/s_nbPeriods;
[~,i]=max(v_meanGamEEG);
s_fdGamEEG=2*i;
v_gamEEGCum= cumsum(v_meanGamEEG);
s_powGamEEG = max(v_gamEEGCum);

NFFT = 2^nextpow2(300/0.36);
f = 3000/2*linspace(0,1,NFFT/2+1);

%% visualization %%
figure(12)
subplot(2,2,1);
    plot(f,v_meanSpVm);
    title('Spectre moyen Spindle vm');
    xlabel('Frequence Hz');
    ylabel('mV^2');
subplot(2,2,2);
    plot(f,v_meanGamVm);
    title('Spectre moyen Gamma vm');
    xlabel('Frequence Hz');
    ylabel('mV^2');
subplot(2,2,3);
    plot(f,v_meanSpEEG);
    title('Spectre moyen spindle EEG');xlabel('Frequence Hz');
    ylabel('mV^2');
subplot(2,2,4);
    plot(f,v_meanGamEEG);
    title('Spectre moyen gamma EEG');xlabel('Frequence Hz');
    ylabel('mV^2');

figure(24)
subplot(2,2,1);
    plot(v_spVmCum);
    title('Spectre cumulatif Spindle vm');xlabel('Frequence Hz');
    ylabel('mV^2');
subplot(2,2,2);
    plot(v_gamVmCum);
    title('Spectre cumulatif Gamma vm');xlabel('Frequence Hz');
    ylabel('mV^2');
subplot(2,2,3);
    plot(v_spEEGCum);
    title('Spectre cumulatif spindle EEG');xlabel('Frequence Hz');
    ylabel('mV^2');
subplot(2,2,4);
    plot(v_gam_EEG_cum);
    title('Spectre cumulatif gamma EEG');xlabel('Frequence Hz');
    ylabel('mV^2');