%computes mean time-frequency map for one cell
%mean is computed by cutting all maps to the same length

prompt = {'Number of periods'};
dlg_title = 'Information required';
num_lines = 1;
def = {'1'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
s_nbPeriods = str2double(answer{1});

m_TFAbsvm1421 = []; %mean of all intra time-frequency maps
m_TFAbsEEG1421 = [];%mean of all EEG time-frequency maps
s_count = 0;

for s_period = 1:s_nbPeriods
    
    str_file = strcat('period', int2str(s_period),'.smr');
    st_info = openSpike2(str_file);
    [m_upstates,v_number,st_info] = upstate_analysis(st_info); 
    
    for s_which = 1:length(v_number)
        
        v_unitVm = m_upstates{s_which,1};
        v_unitEEG = m_upstates{s_which,2};
        b_valid = validity_test(v_unitVm,st_info(1).header.sampleinterval);
        
        if b_valid
            
            s_count = s_count + 1;
            %compute time frequency map
            m_vm = compute_TF(v_unitVm,st_info(1).header.Sampling,20,100,200);
            [m_EEG, v_TimeAxis, v_FreqAxis] = ...
            compute_TF(v_unitEEG,st_info(2).header.Sampling,20,100,200);
            %sum
            if isempty(m_TFAbsvm1421)
                m_TFAbsvm1421 = m_vm;
            else s_W = min(size(m_vm,2),size(m_TFAbsvm1421,2));
                m_TFAbsvm1421 = m_TFAbsvm1421(:,1:s_W) + m_vm(:,1:s_W);
            end
            if isempty(m_TFAbsEEG1421)
                m_TFAbsEEG1421 = m_EEG;
            else s_W = min(size(m_EEG,2),size(m_TFAbsEEG1421,2));
                m_TFAbsEEG1421 = m_TFAbsEEG1421(:,1:s_W) + m_EEG(:,1:s_W);
            
            end

        end
    end
end

m_TFAbsvm1421 = m_TFAbsvm1421/s_count;
m_TFAbsEEG1421 = m_TFAbsEEG1421/s_count;

%% Visualize

figure(4)
subplot(1,2,1)
    f_ImageArray(m_TFAbsvm1421, v_TimeAxis(1:size(m_TFAbsvm1421,2)), v_FreqAxis,'colormap','jet');
    xlabel('Temps (s)');
    ylabel('Fr�quence (Hz)');
    title('Intra')
subplot(1,2,2)
    f_ImageArray(m_TFAbsEEG1421, v_TimeAxis(1:size(m_TFAbsEEG1421,2)), v_FreqAxis,'colormap','jet');
    xlabel('Temps (s)');
    ylabel('Fr�quence (Hz)');
    title('EEG')

