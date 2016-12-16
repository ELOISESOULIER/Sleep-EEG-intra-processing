%Compute mean TF relatively to slow wave phase for evoked and spontaneous
%conditions

st_info = openSpike2('c976.4-fqcy.smr');
b_spont = true;
b_evok = true ;
s_nbBins = 100;
s_binWidth=2*pi/s_nbBins;
s_maxLen=5000;

if b_spont
   v_spont = 1000*st_info(30).data';
end
if b_evok   
   v_evok = 1000*st_info(29).data';
end

time = 0:st_info(2).header.sampleinterval:st_info(2).header.sampleinterval*(length(v_Vm)-1);
m_TFVmSpont = [];
m_TFEEGSpont = [];

s_counts = 0;
if b_spont
    
    for up=1:length(v_spont)
        
        v_stim = find(time>v_spont(up));
        s_stim = v_stim(1);%first point after stimulation
        st_infoSpont = st_info;
        st_infoSpont(1).data = v_Vm(s_stim: min(s_stim+s_maxLen,length(v_Vm)));
        st_infoSpont(2).data = v_EEG(s_stim: min(s_stim+s_maxLen,length(v_EEG)));
        [m_upstates_spont,number,st_infoSpont] =upstate_analysis(st_infoSpont);
        %filter slow wave in intra and EEG
        v_filtSW = f_DesignIIRfilter(st_info(1).header.Sampling,1,1.3,[0.5,200]);
        v_SWVm = f_FilterIIR(st_infoSpont(1).data,v_filtSW);
        v_SWEEG = f_FilterIIR(st_infoSpont(2).data,v_filtSW); 
        v_phaseSWVm = angle(v_SWVm);
        v_unitVmSpont = st_infoSpont(1).data(1:number(2));
        v_unitEEGSpont= st_infoSpont(2).data(1:number(2)); 
        
        v_diff =(v_unitVmSpont(2:end) - v_unitVmSpont(1:end-1))/360*1000;
        b_pa = sum(v_diff >10) <0.5;
        b_len = length(v_unitVmSpont)*0.36>500;
        if b_len && b_pa
            s_counts=s_counts+1;
            [m_vm,~,v_FreqAxis] = f_GaborTransformWait(v_unitVmSpont,st_info(1).header.Sampling,15,100,320);
            m_EEG = f_GaborTransformWait(v_unitEEGSpont,st_info(2).header.Sampling,15,100,320);
            m_reduced_vm_spont = zeros(size(m_vm,1),s_nbBins);
            m_reduced_EEG_spont = zeros(size(m_EEG,1),s_nbBins);
            %on moyenne par rapport � la phase des slow waves
            v_numels=zeros(s_nbBins);
            for index=1:length(v_phaseLoc)-1
                s_bin=1;
                while v_phaseLoc(index) < pi - s_bin*s_binWidth
                        s_bin = s_bin+1;
                end
                m_reduced_vm_spont(:,s_bin)= m_reduced_vm_spont(:,s_bin)+ m_vm(:,index);
                m_reduced_EEG_spont(:,s_bin) = m_reduced_EEG_spont(:,s_bin)+ m_EEG(:,index);
                v_numels(s_bin)=v_numels(s_bin)+1;
            end
            for s_bin=1:s_nbBins
                if v_numels(s_bin)~=0
                    m_reduced_vm_spont(:,s_bin) = m_reduced_vm_spont(:,s_bin)./v_numels(s_bin);
                    m_reduced_EEG_spont(:,s_bin) = m_reduced_EEG_spont(:,s_bin)./v_numels(s_bin);
                else continue
                end
            end 
            if isempty(m_TFVmSpont)
                m_TFVmSpont =m_reduced_vm_spont;
            else m_TFVmSpont = m_TFVmSpont + m_reduced_vm_spont;
            end
             if isempty(m_TFEEGSpont)
                m_TFEEGSpont =m_reduced_EEG_spont;
            else m_TFEEGSpont = m_TFEEGSpont + m_reduced_EEG_spont;
             end

            figure(12)
            subplot(1,3,1);
            f_ImageArray(m_reduced_vm_spont,1:1:s_nbBins, v_FreqAxis,'colormap','jet');
            title('Vm');xlabel('bins de phase');ylabel('Fr�quence');
            subplot(1,3,2);
            f_ImageArray(m_reduced_EEG_spont,1:1:s_nbBins, v_FreqAxis,'colormap','jet');
            title('EEG');xlabel('bins de phase');ylabel('Fr�quence');
            subplot(2,3,3);
            time=0:st_info(1).header.sampleinterval/1000:(length(v_unitVmSpont)-1)*st_info(1).header.sampleinterval/1000;
            plot(time,v_unitVmSpont)
            title('Intra'); xlabel('temps (ms)');ylabel('Fr�quence');
            subplot(2,3,6);
            plot(time,v_unitEEGSpont);
            title('EEG');xlabel('temps (ms)');ylabel('Fr�quence');

            figure(13)
            subplot(1,2,1);
            f_ImageArray(m_TFVmSpont,1:1:s_nbBins, v_FreqAxis,'colormap','jet');
            title('Vm')
            subplot(1,2,2);
            f_ImageArray(m_TFEEGSpont,1:1:s_nbBins, v_FreqAxis,'colormap','jet');
            title('EEG')

       keyboard;
        end
        
    end
end
m_TFVmSpont = m_TFVmSpont/s_counts;      
m_TFEEGSpont = m_TFEEGSpont/s_counts;

figure(1)
subplot(1,2,1)
f_ImageArray(m_TFVmSpont, 1:1:s_nbBins, v_FreqAxis,'colormap','jet');
title('Intra Spontan�');xlabel('bins de phase');ylabel('Fr�quence (Hz)');
subplot(1,2,2)
f_ImageArray(m_TFEEGSpont, 1:1:s_nbBins, v_FreqAxis,'colormap','jet');
title('EEG Spontan�');xlabel('bins de phase');ylabel('Fr�quence (Hz)');


m_TF_vm_evok = [];
m_TF_EEG_evok = [];
s_counte=0;

if b_evok
    for up=1:length(v_evok)
        v_evok  = find(time>v_evok(up));
        info_evok =st_info;
        info_evok(1).data = v_Vm(v_evok(1):min(v_evok(1)+s_maxLen,length(v_Vm)));
        info_evok(2).data = v_EEG(v_evok(1):min(v_evok(1)+s_maxLen,length(v_EEG)));
        [m_upstates_evok,number] =upstate_analysis(info_evok);
        v_unitVmevok = info_evok(1).data(1:number(1)+length(m_upstates_evok{1,1}));
        v_unitEEGevok = info_evok(2).data(1:number(1)+length(m_upstates_evok{1,2}));
        
        v_diff =(v_unitVmevok(2:end) - v_unitVmevok(1:end-1))/360*1000;
        b_pa = sum(v_diff >10) <0.5;
        b_len = length(v_unitVmevok)*0.36>500;
        if b_len && b_pa
            s_counte=s_counte+1;
            m_vm = f_GaborTransformWait(v_unitVmevok,st_info(1).header.Sampling,15,100,320);
            [m_EEG, v_TimeAxis, v_FreqAxis] = ...
            f_GaborTransformWait(v_unitEEG,st_info(2).header.Sampling,15,100,320);
            m_reduced_vm_evok = zeros(size(m_vm,1),s_nbBins);
            m_reduced_EEG_evok = zeros(size(m_EEG,1),s_nbBins);
            %on moyenne par rapport � la phase des slow waves
            v_filtSW =  f_DesignIIRfilter(st_info(1).header.Sampling,2,3);
            v_SWVm = f_FilterIIR(v_unit_EEG,v_filtSW);
            v_analytic_SW = hilbert(v_SWVm);
            v_phaseLoc = angle(v_analytic_SW);
            v_numels=zeros(s_nbBins);
            for index=1:length(v_phaseLoc)-1
                s_bin=1;
                while v_phaseLoc(index) < pi - s_bin*s_binWidth
                    s_bin = s_bin+1;
                end
                m_reduced_vm_evok(:,s_bin)= m_reduced_vm_evok(:,s_bin)+ m_vm(:,index);
                m_reduced_EEG_evok(:,s_bin) = m_reduced_EEG_evok(:,s_bin)+ m_EEG(:,index);
                v_numels(s_bin)=v_numels(s_bin)+1;
            end
            for s_bin=1:s_nbBins
                if v_numels(s_bin)~=0
                    m_reduced_vm_evok(:,s_bin) = m_reduced_vm_evok(:,s_bin)./v_numels(s_bin);
                    m_reduced_EEG_evok(:,s_bin) = m_reduced_EEG_evok(:,s_bin)./v_numels(s_bin);
                else continue
                end
            end 
            if isempty(m_TF_vm_evok)
                m_TF_vm_evok =m_reduced_vm_evok;
            else m_TF_vm_evok = m_TF_vm_evok + m_reduced_vm_evok;
            end
             if isempty(m_TF_EEG_evok)
                m_TF_EEG_evok =m_reduced_EEG_evok;
            else m_TF_EEG_evok = m_TF_EEG_evok + m_reduced_EEG_evok;
             end

            figure(12)
            subplot(1,3,1);
            f_ImageArray(m_reduced_vm_evok,1:1:s_nbBins, v_FreqAxis,'colormap','jet');
            title('Vm');xlabel('bins de phase');ylabel('Fr�quence');
            subplot(1,3,2);
            f_ImageArray(m_reduced_EEG_evok,1:1:s_nbBins, v_FreqAxis,'colormap','jet');
            title('EEG');xlabel('bins de phase');ylabel('Fr�quence');
            subplot(2,3,3);
            time=0:st_info(1).header.sampleinterval/1000:(length(v_unitVmSpont)-1)*st_info(1).header.sampleinterval/1000;
            plot(time,v_unitVmevok)
            title('Intra');xlabel('temps (ms)');ylabel('Fr�quence');
            subplot(2,3,6);
            plot(time,v_unitEEGevok);
            title('EEG');xlabel('temps (ms)');ylabel('Fr�quence');

            figure(13)
            subplot(1,2,1);
            f_ImageArray(m_TF_vm_evok,1:1:s_nbBins, v_FreqAxis,'colormap','jet');
            title('Vm')
            subplot(1,2,2);
            f_ImageArray(m_TF_EEG_evok,1:1:s_nbBins, v_FreqAxis,'colormap','jet');
            title('EEG')
            keyboard;
        end
   end
end


m_TF_vm_evok = m_TF_vm_evok/s_counte;      
m_TF_EEG_evok = m_TF_EEG_evok/s_counte;

figure(2)
subplot(1,2,1)
f_ImageArray(m_TF_vm_evok, 1:1:s_nbBins, v_FreqAxis,'colormap','jet');
title('Intra Evoqu�');xlabel('bins de phase');ylabel('Fr�quence (Hz)');
subplot(1,2,2)
f_ImageArray(m_TF_EEG_evok, 1:1:s_nbBins, v_FreqAxis,'colormap','jet');
title('EEG Evoqu�');xlabel('bins de phase');ylabel('Fr�quence (Hz)');