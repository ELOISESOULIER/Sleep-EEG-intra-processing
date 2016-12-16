%computes mean time frequency map for spontaneous and evoked upstates on n
%cycles preceeding and following stimulation time

%% Preliminary information

%add path to relevant data and toolboxes
addpath(genpath('/home/eloise/Stages/Rythm-ICM/Spike2-data/donnees_stimulees/'));
addpath(genpath('/home/eloise/Stages/Rythm-ICM/Matlab/Toolbox/'));
%name of the files containing info for each condition                                                                                                                                                                                                                                
v_filesEvok = [{'c976.4-fqcy.smr'};{'c1046.5-fqcy.smr'};...
     {'c1054.2-fqcy.smr'};{'c1307-fqcy.smr'};{'c1520.9-fqcy.smr'};...
     {'c2051_fqcy.smr'};{'c1615.5-fqcyevok.smr'};{'c2098.4-fqcyEvok.smr'};...
     {'c1277.6-fqcy2Evok.smr'};{'c1036.2_fqcy.smr'};{'c1176.3-evokfqcy2.smr'}];
 v_filesSpont = [{'c976.4-fqcy.smr'};{'c1046.5-fqcy.smr'};...
     {'c1054.2-fqcy.smr'};{'c1307-fqcy.smr'};{'c1520.9-fqcy.smr'};...
     {'c2051_fqcy.smr'};{'c1615.5-fqcyspont.smr'};{'c2098.4-fqcySpont.smr'};...
     {'c1277.6-fqcySpont.smr'};{'c1036.2_fqcy.smr'};{'c1176.3-spontfqcy.smr'}];
%reference name for the output 
ref = [{'c976'};{'1046'};{'1054'};{'1307'};{'1520'};{'2051'};{'1615'};...
    {'2098'};{'1277'};{'1036'};{'1176'}];

assert(length(v_filesEvok) == length(v_filesSpont));
s_nbFiles = length(v_filesEvok);

params = getParams();                               
s_baseline = params.s_baseline;
s_trialLen = params.s_trialLen;
s_TFreso = params.s_TFreso;
s_nbBins = params.s_nbBins;
s_freqMin = params.s_freqMin;


for s_file = s_nbFiles:s_nbFiles
    
%% compute all upstates for each condition    
    %check if spontaneous and evoked upstates are in the same file or not
    %and subsequently compute upstates for each condition
    b_divided = sum(v_filesEvok{s_file}~=v_filesSpont{s_file}(1:length(v_filesEvok{s_file})))>1;
    if b_divided
        
        st_infoSpont = openSpike2(v_filesSpont{s_file});
        [st_infoSpont,v_VmSpont,v_EEGSpont] = open_subsamp_VmEEG(st_infoSpont);
        v_spont = 1000*st_infoSpont(30).data';
        
        st_infoEvok = openSpike2(v_filesEvok{s_file});
        [st_infoEvok,v_VmEvok,v_EEGEvok] = open_subsamp_VmEEG(st_infoEvok);
        v_evok = 1000*st_infoEvok(29).data';
        
        s_sampleinterval = st_infoSpont(1).header.sampleinterval;
        s_Sampling = st_infoSpont(1).header.Sampling;
        
    else
        
        st_info = openSpike2(v_filesSpont{s_file});
        [st_info,v_Vm,v_EEG] = open_subsamp_VmEEG(st_info);
        v_spont = 1000*st_info(30).data';
        v_evok = 1000*st_info(29).data';
        
        s_sampleinterval = st_info(1).header.sampleinterval;
        s_Sampling = st_info(1).header.Sampling;

    end
   
    
%% compute sum of all the previous upstates for each condition
    m_sumTFSpont = [];
    s_countSpont = 0;
    
    for index = 1:length(v_spont)
        
        if b_divided
             v_time = 0:s_sampleinterval/1000:s_sampleinterval/1000*(length(v_EEGSpont)-1);
        else
             v_time = 0:s_sampleinterval/1000:s_sampleinterval/1000*(length(v_EEG)-1);
        end
        
        v_stim = find(v_time>v_spont(index)); 
        s_stim = v_stim(1) - floor(s_baseline/s_sampleinterval*1000);% first point after stimulation+baseline

   
        %check if there are 4 cycles before next stimulation
        if (index == length(v_spont) && (v_time(end)-v_spont(end)) > 2000) || ...
             (index < length(v_spont) &&  index > 1 &&...
             (v_spont(index+1) - v_spont(index))> 2000 ...
             && (v_spont(index) - v_spont(index-1)) > 2000)  || ...
             (index == 1 && (v_spont(index+1) - v_spont(index))>2000)
            
            %count number of samples taken into account in the mean
            s_countSpont = s_countSpont + 1;
            
            if b_divided
                v_signal = v_EEGSpont(s_stim:s_stim+floor(s_trialLen/s_sampleinterval*1000));
                [m_TF,v_TimeAxis,v_FreqAxis] = compute_TF(v_signal,s_Sampling,...
                    params.s_freqMin,100,params.s_TFreso);
            else
                 v_signal = v_EEG(s_stim:s_stim+floor(s_trialLen/s_sampleinterval*1000));
                [m_TF,v_TimeAxis,v_FreqAxis] = compute_TF(v_signal,s_Sampling,...
                    params.s_freqMin,100,params.s_TFreso);
            end

            if isempty(m_sumTFSpont)
                m_sumTFSpont = m_TF;
            else
                m_sumTFSpont = m_sumTFSpont + m_TF;
            end


        else continue

        end


    end

    m_sumTFSpont = m_sumTFSpont/s_countSpont;

    m_sumTFEvok = [];
    s_countEvok = 0;
    
    for index = 1:length(v_evok)

        if b_divided
             v_time = 0:s_sampleinterval/1000:s_sampleinterval/1000*(length(v_EEGEvok)-1);
        else
             v_time = 0:s_sampleinterval/1000:s_sampleinterval/1000*(length(v_EEG)-1);
        end
        
        v_stim = find(v_time>v_evok(index)); 
        s_stim = v_stim(1) - floor(s_baseline/s_sampleinterval*1000);% first point after stimulation+baseline

        %if there are 4 cycles before next stimulation
        if (index == length(v_evok) && (v_time(end)-v_evok(end))>2000) || ...
             (index < length(v_evok) &&  index > 1 &&...
             (v_evok(index+1) - v_evok(index))>2000 ...
             && (v_evok(index) - v_evok(index-1)) > 2000)  || ...
             (index == 1 && (v_evok(index+1) - v_evok(index))>2000)
            
            s_countEvok = s_countEvok + 1;
            
            if b_divided
                v_signal = v_EEGEvok(s_stim:s_stim+floor(s_trialLen/s_sampleinterval*1000));
                [m_TF,v_TimeAxis,v_FreqAxis] = compute_TF(v_signal,s_Sampling,...
                    params.s_freqMin,100,params.s_TFreso);
            else
                 v_signal = v_EEG(s_stim:s_stim+floor(s_trialLen/s_sampleinterval*1000));
                [m_TF,v_TimeAxis,v_FreqAxis] = compute_TF(v_signal,s_Sampling,...
                    params.s_freqMin,100,params.s_TFreso);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 params.s_freqMin,100,params.s_TFreso);
            end

            if isempty(m_sumTFEvok)
                m_sumTFEvok = m_TF;
            else
                m_sumTFEvok = m_sumTFEvok + m_TF;
            end


        else continue

        end


    end

    m_sumTFEvok = m_sumTFEvok/s_countEvok;

%% Visualisation    
    
    figure(1)
    subplot(1,2,1)
        v_clim = [prctile(m_sumTFSpont(:),1),prctile(m_sumTFSpont(:),100)];
        f_ImageArray(m_sumTFSpont,v_TimeAxis,v_FreqAxis,'colormap','jet','limits',v_clim)
        title(['Spontané',int2str(s_countSpont)])
        colorbar
    subplot(1,2,2)
        f_ImageArray(m_sumTFEvok,v_TimeAxis,v_FreqAxis,'colormap','jet','limits',v_clim)
        title(['Evoqué',int2str(s_countEvok)])
        colorbar
        saveas(gcf,['/home/eloise/Stages/Rythm-ICM/Matlab/cycle/',ref{s_file},'.png']);

%% save output under chosen ref name 
    save(['/home/eloise/Stages/Rythm-ICM/Matlab/cycle/',ref{s_file}, '-Spont'],'m_sumTFSpont');
    save(['/home/eloise/Stages/Rythm-ICM/Matlab/cycle/',ref{s_file} ,'-Evok'],'m_sumTFEvok');

end