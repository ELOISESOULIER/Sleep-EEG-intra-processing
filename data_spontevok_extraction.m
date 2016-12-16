% Extracting data for spontaneous and evoked conditions in .mat
%exactly like data_extraction but for spontaneous and evoked (could be
%optimized...)

files_evok = [{'c976.4-fqcy.smr'};{'c1036.2_fqcy.smr'};{'c1046.5-fqcy.smr'};...
 {'c1054.2-fqcy.smr'};{'c1307-fqcy.smr'};{'c1520.9-fqcy.smr'};...
 {'c2051_fqcy.smr'};{'c1176.3-evokfqcy1.smr'};{'c1277.6-fqcy2evok.smr'};...
 {'c1615.5-fqcyevok.smr'};{'c2098.4-fqcyEvok.smr'}];
files_spont = [{'c976.4-fqcy.smr'};{'c1036.2_fqcy.smr'};{'c1046.5-fqcy.smr'};...
 {'c1054.2-fqcy.smr'};{'c1307-fqcy.smr'};{'c1520.9-fqcy.smr'};...
 {'c2051_fqcy.smr'};{'c1176.3-spontfqcy.smr'};{'c1277.6-fqcySpont.smr'};...
 {'c1615.5-fqcyspont.smr'};{'c2098.4-fqcySpont.smr'}];

assert(length(files_spont) == length(files_evok));
nb_bins=60;
s_nbFiles = length(files_evok); 

tfspont = cell(1,2);
tfevok = cell(1,2);
s_countspont=0;
s_countevok = 0;

for file = 1:s_nbFiles
    
%% open files and extract info
    info = openSpike2(files_spont{file});
    [m_spont,~,~] = spont_evok_extract(info,true,false);
    info = openSpike2(files_evok{file});
    [~,m_evok,info] = spont_evok_extract(info,false,true);
    Si =info(1).header.sampleinterval;
    
%% Compute TF map 
    
    %spontaneous
    for s_up = 1:length(m_spont)
        
        b_valid = validity_test(m_spont{s_up,1},Si);
        
        if b_valid
            
            s_countspont = s_countspont+1;
            v_unitVmSpont=m_spont{s_up,1};
            v_unitEEGSpont=m_spont{s_up,2};
            
            m_vmSpont = compute_TF(v_unitVmSpont,info(1).header.Sampling,...
                20,80,60);
            [m_EEGSpont, v_TimeAxis, v_FreqAxis] = compute_TF(...
                v_unitEEGSpont,info(2).header.Sampling,20,80,60);
            
            %compute reduced map by averaging over given number of bins
            m_redVmSpont = zeros(size(m_vmSpont,1),nb_bins);
            m_redEEGSpont = zeros(size(m_EEGSpont,1),nb_bins);
            s_W = size(m_vmSpont,2);
            s_WBin =floor(s_W/nb_bins);
            v_aux = 1:s_WBin:nb_bins*s_WBin;
            v_TimeAxis = v_TimeAxis(v_aux);
            for i=1:nb_bins
                m_redEEGSpont(:,i) = mean(m_EEGSpont(:,(i-1)*s_WBin +1: i*s_WBin),2);
                m_redVmSpont(:,i) = mean(m_vmSpont(:,(i-1)*s_WBin +1: i*s_WBin),2);
            end
            tfspont{s_countspont,1}=m_redVmSpont;
            tfspont{s_countspont,2}=m_redEEGSpont;
            
        end
        
    end
    
    for s_up = 1:length(m_evok)
       
        b_valid = validity_test(m_evok{s_up,1});
        
        if b_valid
            
            s_countevok=s_countevok+1;
            v_unitVmEvok = m_evok{s_up,1};
            v_unitEEGEvok = m_evok{s_up,2};
            
            m_vmEvok = compute_TF(v_unitVmEvok,info(1).header.Sampling,20,80,60);
            [m_EEGEvok, v_TimeAxis, v_FreqAxis] = ...
            compute_TF(v_unitEEGEvok,info(2).header.Sampling,20,80,60);
        
            m_redVmevok = zeros(size(m_vmEvok,1),nb_bins);
            m_redEEGevok = zeros(size(m_EEGEvok,1),nb_bins);
            s_W = size(m_vmEvok,2);
            s_WBin =floor(s_W/nb_bins);
            v_aux = 1:s_WBin:nb_bins*s_WBin;
            v_TimeAxis = v_TimeAxis(v_aux);
            for i=1:nb_bins
                m_redEEGevok(:,i) = mean(m_EEGEvok(:,(i-1)*s_WBin +1: i*s_WBin),2);
                m_redVmevok(:,i) = mean(m_vmEvok(:,(i-1)*s_WBin +1: i*s_WBin),2);
            end
            tfevok{s_countevok,1}=m_redVmevok;
            tfevok{s_countevok,2}=m_redEEGevok;
            
        end
        
    end
    
end

save('tfspont','tfspont');
save('tfevok','tfevok');
