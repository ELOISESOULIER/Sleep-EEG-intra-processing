%computes power spectrum for spontaneous and evoked upstates for one single
%cell

str_cellNameEvok = 'c2098.4-fqcyEvok.smr';
str_cellNameSpont = 'c2098.4-fqcySpont.smr';

%read data and extract upstates
st_info = openSpike2(str_cellNameEvok);
[~,m_evok,st_info] = spont_evok_extract(st_info,false,true);
st_info =  openSpike2(str_cellNameSpont);
[m_spont,~,st_info] = spont_evok_extract(st_info,true,false);

%% compute mean spectrum over all upstates of the cell

%spontaneous
[v_spVmSpont, v_gamVmSpont, v_spEEGSpont,v_gamEEGSpont]...
    = mean_spec(m_spont);
v_spVmCumSpont= cumsum(v_spVmSpont); s_powSpVmSpont = max(v_spVmCumSpont);
v_gamVmCumSpont= cumsum(v_gamVmSpont); s_powGamVmSpont = max(v_gamVmCumSpont);
v_spEEGCumSpont= cumsum(v_spEEGSpont); s_powSpEEGSpont = max(v_spEEGCumSpont);
v_gamEEGCumSpont = cumsum(v_gamEEGSpont); s_powGamEEGSpont = max(v_gamEEGCumSpont);
[~,i]=max(v_spVmSpont);s_fdSpVmSpont=2*i;
[~,i]=max(v_gamVmSpont);s_fdGamVmSpont=2*i;
[~,i]=max(v_spEEGSpont);s_fdSpEEGSpont=2*i;
[~,i]=max(v_gamEEGSpont);s_fdGamEEGSpont=2*i;  

%evoked
[v_spVmEvok, v_gamVmEvok, v_spEEGEvok,v_gamEEGEvok]...
    = mean_spec(m_evok);
v_spVmCumEvok= cumsum(v_spVmEvok); pow_sp_vm_evok = max(v_spVmCumEvok);
v_gamVmCumEvok= cumsum(v_gamVmEvok); pow_gam_vm_evok = max(v_gamVmCumEvok);
v_spEEGCumEvok= cumsum(v_spEEGEvok); pow_sp_eeg_evok = max(v_spEEGCumEvok);
v_gamEEGCumEvok = cumsum(v_gamEEGEvok); pow_gam_EEG_evok = max(v_gamEEGCumEvok);
[~,i]=max(v_spVmEvok);s_fdSpVmEvok=2*i;
[~,i]=max(v_gamVmEvok);s_fdGamVmEvok=2*i;
[~,i]=max(v_spEEGEvok);s_fdSpEEGEvok=2*i;
[~,i]=max(v_gamEEGEvok);s_fdGamEEGEvok=2*i;

%% visualisation

NFFT = 2^nextpow2(300/0.36);
f = 3000/2*linspace(0,1,NFFT/2+1);

figure(12)
    subplot(2,2,1);
    plot(f,v_spVmSpont,'b');
    hold on;
    plot(f,v_spVmEvok,'r');
%     hold off;
    legend('Spontan�','Evoqu�');
    title('Spectre moyen Spindle vm');
    xlabel('Frequence Hz');
    ylabel('mV^2');
subplot(2,2,2);
    plot(f,v_gamVmSpont,'b');
    hold on;
    plot(f,v_gamVmEvok,'r');
%     hold off;
    legend('Spontan�','Evoqu�');
    title('Spectre moyen Gamma vm');
    xlabel('Frequence Hz');
    ylabel('mV^2');
subplot(2,2,3);
    plot(f,v_spEEGSpont,'b');
    hold on;
    plot(f,v_spEEGEvok,'r');
%     hold off;
    legend('Spontan�','Evoqu�');
    title('Spectre moyen spindle EEG ');xlabel('Frequence Hz');
    ylabel('muV^2');
subplot(2,2,4);
    plot(f,v_gamEEGSpont,'b');
    hold on;
    plot(f,v_gamEEGEvok,'r');
%     hold off;
    legend('Spontan�','Evoqu�');
    title('Spectre moyen gamma EEG ');xlabel('Frequence Hz');
    ylabel('muV^2');
