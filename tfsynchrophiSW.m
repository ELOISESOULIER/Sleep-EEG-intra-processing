%compute mean phase-frequency map for one cell : 
%mean is computed by reducing all maps to a given number of bins based on
%the value of  the corresponding slow wave angle

s_nbPeriods = 1;
s_nbBins = 200;
s_binWidth = 2*pi/s_nbBins;
v_phaseAxis = -pi:s_binWidth:-pi+s_binWidth*(s_nbBins-1);
params = getParams();

m_phiFVm1046 = []; %mean phase-frequency map for intra
m_phiFEEG1046 = []; %mean phase-frequency map for EEG
v_meanEEG1046 = zeros(1,s_nbBins);
v_meanVm1046 = zeros(1,s_nbBins);
v_numelsTotVm = zeros(1,s_nbBins); 
v_numelsTotEEG = zeros(1,s_nbBins);
s_count = 0;

for s_period = 1:s_nbPeriods
    
    %read files and extract upstates
    str_file = strcat('period', int2str(s_period),'.smr');
    st_info = openSpike2(str_file);
    [m_upstates,v_number,st_info] = upstate_analysis(st_info); 
    v_Vm = st_info(1).data;
    v_EEG  = st_info(2).data;
    s_sampleInt = st_info(1).header.sampleinterval;
    
    %compute slow waves
    v_filtSW =  f_DesignIIRfilter(st_info(1).header.Sampling,2,3);
    v_SWEEG = f_FilterIIR(v_EEG,v_filtSW);
    v_SWVm = f_FilterIIR(v_Vm,v_filtSW);
    v_phaseSWEEG = angle( hilbert(v_SWEEG) );
    v_phaseSWVm = angle( hilbert(v_SWVm - mean(v_SWVm)) );

    for s_which = 1:length(v_number)-1  
%% compute time-frequency map 
        v_unitVm = v_Vm(max(1,v_number(s_which)):v_number(s_which+1));
        v_unitEEG = v_EEG(max(1,v_number(s_which)):v_number(s_which+1));
        v_phaseLocVm = v_phaseSWVm(max(1,v_number(s_which)):v_number(s_which+1));
        v_phaseLocEEG = v_phaseSWEEG(max(1,v_number(s_which)):v_number(s_which+1));
        b_valid = validity_test(v_unitVm, st_info(1).header.sampleinterval);
        if b_valid
            s_count = s_count + 1;

            %time-frequency map
            m_vm = compute_TF(v_unitVm,st_info(1).header.Sampling,...
                params.s_freqMin,params.s_freqMax,params.s_TFreso);
            [m_EEG,~, v_FreqAxis] = compute_TF(v_unitEEG,...
                st_info(2).header.Sampling,params.s_freqMin,...
                params.s_freqMax,params.s_TFreso);

 %% compute phase-frequency map
            %reduce to s_nbBins bins according to phase
            m_redVm = zeros(size(m_vm,1),s_nbBins);
            m_redEEG = zeros(size(m_EEG,1),s_nbBins);
            v_redVm = zeros(1,s_nbBins);
            v_redEEG = zeros(1,s_nbBins);
            v_numelsVm = zeros(1,s_nbBins);
            v_numelsEEG = zeros(1,s_nbBins);
            %sort by slow-wave angle value
            for index = 1:length(v_phaseLocVm)-1
                s_bin = 1;
                while v_phaseLocVm(index) > -pi + s_bin*s_binWidth
                    s_bin = s_bin+1;
                end
                m_redVm(:,s_bin) = m_redVm(:,s_bin)+ m_vm(:,index);
                v_redVm(s_bin) = v_redVm(s_bin) + v_unitVm(index);
                v_numelsVm(s_bin) = v_numelsVm(s_bin)+1;
            end
            for index = 1:length(v_phaseLocEEG)-1
                s_bin = 1;
                while v_phaseLocEEG(index) > -pi + s_bin*s_binWidth
                    s_bin = s_bin+1;
                end
                m_redEEG(:,s_bin) = m_redEEG(:,s_bin)+ m_EEG(:,index);
                v_redEEG(s_bin) = v_redEEG(s_bin) + v_unitEEG(index);
                v_numelsEEG(s_bin) = v_numelsEEG(s_bin)+1;
            end
            %compute average according to number of elements in each bin
            for s_bin = 1:s_nbBins
                if v_numelsVm(s_bin) ~= 0
                    m_redVm(:,s_bin) = m_redVm(:,s_bin)./v_numelsVm(s_bin);
                    v_redVm(s_bin) = v_redVm(s_bin)/v_numelsVm(s_bin);
                    v_numelsVm(s_bin) = 1;
                else continue
                end
            end 
            for s_bin = 1:s_nbBins
                if v_numelsEEG(s_bin) ~= 0
                    m_redEEG(:,s_bin) = m_redEEG(:,s_bin)./v_numelsEEG(s_bin);
                    v_redEEG(s_bin) = v_redEEG(s_bin)/v_numelsEEG(s_bin);
                    v_numelsEEG(s_bin) = 1;
                else continue
                end
            end 
            %compute sum of all maps
            for s_bin=1:s_nbBins
                if isempty(m_phiFVm1046)
                    m_phiFVm1046 = m_redVm;
                else m_phiFVm1046 = m_phiFVm1046 + m_redVm;
                end
            end
            v_meanVm1046 = v_meanVm1046 + v_redVm;
             if isempty(m_phiFEEG1046)
                m_phiFEEG1046 = m_redEEG;
            else m_phiFEEG1046 = m_phiFEEG1046 + m_redEEG;
             end
            v_meanEEG1046 = v_meanEEG1046 + v_redEEG;
            
            v_numelsTotVm = v_numelsTotVm + v_numelsVm;
            v_numelsTotEEG = v_numelsTotEEG + v_numelsEEG;
            
            %visualisation
%             figure(12)
%             time=0:s_sampleInt/1000:(length(v_unitVm)-1)*s_sampleInt/1000;
%             subplot(2,2,3);
%                 f_ImageArray(m_redVm,v_phaseAxis, v_FreqAxis,'colormap','jet');
%                 title('Vm');xlabel('Phase des ondes lentes');ylabel('Fr�quence');
%             subplot(2,2,4);
%                 f_ImageArray(m_redEEG,v_phaseAxis, v_FreqAxis,'colormap','jet');
%                 title('EEG');xlabel('Phase des ondes lentes');ylabel('Fr�quence');
%             subplot(2,2,1);
%                 plot(v_phaseAxis,v_redVm);
%                 title('Intra');xlabel('Phase (rad)');ylabel('Fr�quence');
%             subplot(2,2,2);
%                 plot(v_phaseAxis,v_redEEG);
%                 title('EEG');xlabel('Phase (rad)');ylabel('Fr�quence');
%             
%             figure(13)
%             subplot(2,2,3);
%                 f_ImageArray(m_TFVm,v_phaseAxis, v_FreqAxis,'colormap','jet');
%                 title('Vm'); xlabel('Phase des ondes lentes'); ylabel('Fr�quence');
%             subplot(2,2,4);
%                 f_ImageArray(m_TFEEG,v_phaseAxis, v_FreqAxis,'colormap','jet');
%                 title('EEG'); xlabel('Phase des ondes lentes'); ylabel('Fr�quence');
%             subplot(2,2,1)
%                  plot(v_phaseAxis,v_meanVm);
%                 title('Intra');xlabel('Phase (rad)');ylabel('Fr�quence');
%             subplot(2,2,2);
%                 plot(v_phaseAxis,v_meanEEG);
%                 title('EEG');xlabel('Phase (rad)');ylabel('Fr�quence');
%        keyboard;

        end
    
    end
end

%% Compute average of all maps
for s_bin = 1:s_nbBins
    if v_numelsTotVm(s_bin) ~= 0
        m_phiFVm1046(:,s_bin) = m_phiFVm1046(:,s_bin)./v_numelsTotVm(s_bin);
    else continue
    end
end 
v_meanVm1046 = v_meanVm1046./v_numelsTotVm;
for s_bin=1:s_nbBins
    if v_numelsTotEEG(s_bin) ~= 0
        m_phiFEEG1046(:,s_bin) = m_phiFEEG1046(:,s_bin)./v_numelsTotEEG(s_bin);
    else continue
    end
end
v_meanEEG1046 = v_meanEEG1046./v_numelsTotVm;

%% Visualize result %%
figure(1)
subplot(2,2,1)
    f_ImageArray(m_phiFVm1046,v_phaseAxis, v_FreqAxis,'colormap','jet');
    title('Intra');xlabel('Phase des ondes lentes (rad)');ylabel('Fr�quence (Hz)');
subplot(2,2,2)
    f_ImageArray(m_phiFEEG1046, v_phaseAxis, v_FreqAxis,'colormap','jet');
    title('EEG');xlabel('Phase des ondes lentes (rad)');ylabel('Fr�quence (Hz)');
subplot(3,2,5)
    plot(v_phaseAxis,v_meanVm1046);
    title('Intra');xlabel('Phase (rad)');ylabel('Potentiel(mV)');
subplot(3,2,6);
    plot(v_phaseAxis,v_meanEEG1046);
    title('EEG');xlabel('Phase (rad)');ylabel('Potentiel(muV)');


