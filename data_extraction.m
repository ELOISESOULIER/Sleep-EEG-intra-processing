% Extract time frequency maps for each upstate in .mat so as to use them in
% python

s_nbPeriods=6;
tf_1036 = cell(1,2);
s_num = 1036;
s_nbBins=20;

s_count=0;

for s_period=1:s_nbPeriods
    
    %open file and extract upstates
    str_file = strcat('period', int2str(s_period),'.smr');
    st_info = openSpike2(str_file);
    [m_upstates,number,st_info]=upstate_analysis(st_info);
    
    %for each upstate, compute reduced map and store it
    for s_up = 1:length(m_upstates)
        
        v_unitVm=m_upstates{s_up,1};
        v_unitEEG=m_upstates{s_up,2};
        b_valid = validity_test(v_unitVm,st_info(1).header.sampleinterval);
        
        if b_len && b_pa
            
            s_count = s_count + 1;   
            [m_vm, v_TimeAxis, v_FreqAxis] = compute_TF(v_unitVm,...
                st_info(1).header.Sampling,20,80,60);
            m_EEG =compute_TF(v_unitEEG,st_info(2).header.Sampling,20,80,60);
            
            %reducing map to given number of bins by averaging
            m_redVm = zeros(size(m_vm,1),s_nbBins);
            m_redEEG = zeros(size(m_EEG,1),s_nbBins);
            s_W = size(m_vm,2);
            s_WBin =floor(s_W/s_nbBins);
            v_aux = 1:s_WBin:s_nbBins*s_WBin;
            v_TimeAxis = v_TimeAxis(v_aux);
            for i=1:s_nbBins
                m_redEEG(:,i) = mean(m_EEG(:,(i-1)*s_WBin +1: i*s_WBin),2);
                m_redVm(:,i) = mean(m_vm(:,(i-1)*s_WBin +1: i*s_WBin),2);
            end
            tf_1036{s_count,1}=m_redVm;
            tf_1036{s_count,2}=m_redEEG;
            
        end
        
    end
    
end

save(['tf_',int2str(s_num)],['tf_',int2str(s_num)]);