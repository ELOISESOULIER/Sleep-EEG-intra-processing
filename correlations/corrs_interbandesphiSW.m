nb_period=5;
freq_min=10;
freq_max=60;
corr_bin=8;
nb_bins=100;
s_binWidth=2*pi/nb_bins;
nb_bands=100;%2*(freq_max-freq_min):


for s_period =1:nb_period
    
    file = strcat('period', int2str(s_period),'.smr');
    st_info = openSpike2(file);
    [m_upstates,number,st_info] = upstate_analysis(st_info); 
    v_vm=st_info(1).data;
    v_EEG = st_info(2).data;
    %récupérer l'onde lente
    v_filtSW =  f_DesignIIRfilter(st_info(1).header.Sampling,2,3);
    v_SW = f_FilterIIR(v_EEG,v_filtSW);
    v_analytic_SW = v_SW+hilbert(v_SW)*i;
    v_phaseSW = angle(v_analytic_SW);

    %découpe les upstates dans EEG et intra

    s_which=1;
    while ~isempty(m_upstates{s_which,1}) && s_which<length(number)-1
        v_unitVm = v_vm(number(s_which):number(s_which+1));
        v_unitEEG =  v_EEG(number(s_which):number(s_which+1));
        %on vérifie longueur et absence de pa
        b_valid = validity_test(v_unitVm,st_info(1).header.sampleinterval);
        if b_valid
            %on établit un tf pour vm et EEG
            [m_vm, v_TimeAxis, v_FreqAxis] = ...
            f_GaborTransformWait(v_unitVm,st_info(1).header.Sampling,freq_min,freq_max,nb_bands);
            [m_EEG, v_TimeAxis, v_FreqAxis] = ...
            f_GaborTransformWait(v_unitEEG,st_info(2).header.Sampling,freq_min,freq_max,nb_bands);
            
            %on moyenne par rapport à la phase des slow waves
            m_reduced_vm = zeros(size(m_vm,1),nb_bins);
            m_reduced_EEG = zeros(size(m_EEG,1),nb_bins);
            v_phaseLoc = v_phaseSW(number(s_which):number(s_which)+length(v_unitVm)-1);
            v_numels=zeros(1,nb_bins);
            for index=1:length(v_phaseLoc)-1
                s_bin=1;
                while v_phaseLoc(index) < pi - s_bin*s_binWidth
                    s_bin = s_bin+1;
                end
                m_reduced_vm(:,s_bin)= m_reduced_vm(:,s_bin)+ m_vm(:,index);
                m_reduced_EEG(:,s_bin) = m_reduced_EEG(:,s_bin)+ m_EEG(:,index);
                v_numels(s_bin)=v_numels(s_bin)+1;
            end
            v_numels(v_numels==0)=1;
            for s_bin=1:nb_bins
                m_reduced_vm(:,s_bin) = m_reduced_vm(:,s_bin)./v_numels(s_bin);
                m_reduced_EEG(:,s_bin) = m_reduced_EEG(:,s_bin)./v_numels(s_bin);
            end 
            
            %on construit une matrice des corrélations interbandes
            M_aux=zeros(nb_bands);
            for abs=1:nb_bands
                v_vmref=m_reduced_vm(abs,:) - mean(m_reduced_vm(abs,:));
                for ord=1:nb_bands
                    M_aux(abs,ord)=max(xcorr(v_vmref,m_reduced_EEG(ord,:)-mean(m_reduced_EEG(ord,:)),'coeff'));
                end
            end
            M=zeros(floor(nb_bands/corr_bin));
            for abs=1:floor(nb_bands/corr_bin)
                for ord=1:floor(nb_bands/corr_bin)
                    coeff=[];
                    for index=1:corr_bin
                        coeff= [coeff,M_aux((abs-1)*corr_bin+index,(ord-1)*corr_bin+index)];
                    M(abs,ord)=mean(coeff);
                    end
                end
            end
            
            pas = (freq_max - freq_min)/nb_bands*corr_bin;
            v=freq_min+pas:pas:freq_max;
            
            s_W = size(m_vm,2);
            s_WBin =floor(s_W/nb_bins);
            v_aux = 1:s_WBin:nb_bins*s_WBin;
            v_TimeAxis = v_TimeAxis(v_aux);
           
            figure(7)
            subplot(1,2,1)
            f_ImageArray(m_reduced_vm,v_TimeAxis, v_FreqAxis,'colormap','jet');
            title('Vm');
            subplot(1,2,2);
            f_ImageArray(m_reduced_EEG,v_TimeAxis, v_FreqAxis,'colormap','jet')
            title('EEG');
            
            figure(8)
            title('Correlations');
            f_ImageArray(M,v,v,'colormap','jet');
            
            keyboard;
        end
        s_which=s_which+1;
    end
end
                    
                    