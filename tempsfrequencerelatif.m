%compute mean time-frequency map for one cell
%mean is computed by reducing all maps to a given number of bins

prompt = {'Number of periods'};
dlg_title = 'Information required';
num_lines = 1;
def = {'1'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
s_nbPeriods = str2num(answer{1});

params = getParams();
m_TFRelvm1421 = []; %sum of all intra time-frequency maps
m_TFRelEEG1421 = []; %sum of all EEG time-frequency maps
s_nbBins = 400;
s_count = 0;

for s_period =1:s_nbPeriods
    
    st_file = strcat('period', int2str(s_period),'.smr');
    str_info = openSpike2(st_file);
    [m_upstates,v_number,str_info] = upstate_analysis(str_info); 
        
    for s_which = 1:length(number)
        
        v_unitVm =m_upstates{s_which,1};
        v_unitEEG = m_upstates{s_which,2};
        b_valid = validity_test(v_unitVm,str_info(1).header.sampleinterval);
        
        if b_valid
            
            s_count = s_count + 1;
       
            %time-frequency map
            m_vm = compute_TF(v_unitVm,str_info(1).header.Sampling,...
                params.s_freqMin,params.s_freqMax,params.s_TFreso);
            [m_EEG, v_TimeAxis, v_FreqAxis] = ...
            compute_TF(v_unitEEG,str_info(2).header.Sampling,...
                params.s_freqMin,params.s_freqMax,params.s_TFreso);

            %reduce to s_nbBins bins
            m_redVm = zeros(size(m_vm,1),s_nbBins);
            m_redEEG = zeros(size(m_EEG,1),s_nbBins);
            s_WBin =floor(size(m_vm,2)/s_nbBins);
            v_aux = 1:s_WBin:s_nbBins*s_WBin;
            v_TimeAxis = v_TimeAxis(v_aux);
            for i=1:s_nbBins
                m_redEEG(:,i) = mean(m_EEG(:,(i-1)*s_WBin +1: i*s_WBin),2);
                m_redVm(:,i) = mean(m_vm(:,(i-1)*s_WBin +1: i*s_WBin),2);
            end
            
            % compute sum
            if isempty(m_TFRelvm1421)
                m_TFRelvm1421 =m_redVm;
            else m_TFRelvm1421 = m_TFRelvm1421 + m_redVm;
            end
            if isempty(m_TFRelEEG1421)
                m_TFRelEEG1421 =m_redEEG;
            else m_TFRelEEG1421 = m_TFRelEEG1421 + m_redEEG;
            end

        end
    end
end

m_TFRelvm1421 = m_TFRelvm1421/s_count;
m_TFRelEEG1421 = m_TFRelEEG1421/s_count;

%% Visualize

figure(1)
subplot(1,2,1)
    f_ImageArray(m_TFRelvm1421, v_TimeAxis, v_FreqAxis,'colormap','jet');
    xlabel('Temps (%upstate)');
    ylabel('Fr�quence (Hz)');
    title('Intra')  
subplot(1,2,2)
    f_ImageArray(m_TFRelEEG1421, v_TimeAxis, v_FreqAxis,'colormap','jet');%,'limits',v_Clim)
    xlabel('Temps (%upstate)');
    ylabel('Fr�quence (Hz)');
    title('EEG')  

