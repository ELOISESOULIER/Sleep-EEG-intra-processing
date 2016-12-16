%extracts spontaneous and evoked upstates from Spike2 file containing intra
%EEG and stimulation times
%st_info : file opened with openSpike2 containine intra as channel 1, EEG
%as channel 2, and nonstim times as channel 30 or stim times as channel 29
%b_spont (resp b_evok) : boolean specifying whether spontaneous (resp 
%evoked) upstates exist
%b_evok : boolean specifying whether evoked upstates ae wanted as output
%b_total : default value = false, only the upstate signal is returned
% if true, 0.1 sec before stimulation time and 1 sec after stimulation time
% is returned
%m_spont (resp m_evok) : 2-column cell containing the spontaneous (resp 
%evoked) upstates and the corresponding parts of the EEG

function [m_spont,m_evok,st_newInfo] = spont_evok_extract(st_info,b_spont,b_evok,b_total)
 
    params = getParams();
    
    if nargin == 1
        b_spont = true;
        b_evok =true;
    end
    
    if nargin == 1 || nargin == 3
        b_total = false;
    end
    
    if nargin == 0 || nargin == 2
        error('ERROR: Invalid number of arguments')
    end

    if b_spont
        v_spont = 1000*st_info(30).data'; %to ms
        s_numUps = length(v_spont);
    else m_spont=0;
    end
    
    if b_evok   
        v_evok = 1000*st_info(29).data';
        s_numUps = length(v_evok);
    else m_evok = 0;
    end

    s_maxLen = 6000;
    s_baseline = params.s_baseline;
    s_trialLen = params.s_trialLen;
    v_VmOrigin = st_info(1).data;
    v_EEG = st_info(2).data;
    %setting intra and EEG to the same sampling rate
    s_subsampRate = st_info(2).header.sampleinterval/st_info(1).header.sampleinterval;
    v_index = 1:s_subsampRate:length(v_VmOrigin);
    st_info(1).header.Sampling = st_info(1).header.Sampling/s_subsampRate;
    st_info(1).header.sampleinterval = st_info(1).header.sampleinterval*s_subsampRate;
    st_newInfo = st_info;
    v_Vm = v_VmOrigin(v_index);
    %assert equal length
    s_min = min(length(v_EEG),length(v_Vm));
    v_Vm = v_Vm(1:s_min); 
    v_EEG = v_EEG(1:s_min);
    st_info(1).data = v_Vm;
    st_info(2).data = v_EEG;
    assert(length(v_Vm)==length(v_EEG));
%% Visualize position of stimulations on whole signal    

    v_time = 0:st_info(1).header.sampleinterval/1000:st_info(1).header.sampleinterval/1000*(length(v_Vm)-1);

%     figure()
%     plot(v_time,v_Vm)
%     hold on;
%     if b_spont
%         plot(v_spont,-60*ones(1,s_numUps),'r*')
%          hold on;
%     end
%     if b_evok
%         plot(v_evok,-60*ones(1,s_numUps),'y*');
%     end
%     hold off;
                                                                                                                        
%% spontaneous %%    
    if b_spont
        m_spont = cell(s_numUps,2);
        for index = 1:s_numUps
            v_stim = find(v_time>v_spont(index)); 
            s_stim = v_stim(1);% first point after stimulation
            if b_total
                %add baseline before stim
                s_stim = s_stim - floor(s_baseline/(st_info(1).header.sampleinterval/1000));
            end
            st_infoSpont = st_info;
            st_infoSpont(1).data = v_Vm(s_stim: min(s_stim+s_maxLen,end));
            st_infoSpont(2).data = v_EEG(s_stim: min(s_stim+s_maxLen,end));
            [m_upstatesSpont,v_number] = upstate_analysis(st_infoSpont);
            if ~b_total %cut on the edges of the upstate
                m_spont{index,1} = m_upstatesSpont{1,1};
                m_spont{index,2} = m_upstatesSpont{1,2};
            else  m_spont{index,1} = st_infoSpont(1).data(1:...
                    min(floor((s_baseline + s_trialLen)/(st_info(1).header.sampleinterval/1000)),end));
                  m_spont{index,2} = st_infoSpont(2).data(1:...
                      min(floor((s_baseline + s_trialLen)/(st_info(1).header.sampleinterval/1000)),end));
                  
            end
            %check that "upstate_analysis" did not select the second
            %upstate after stim time
            if ~b_total
                if v_number(1)>length(st_infoSpont(1).data)/10
                    m_spont{index,1} = 0;
                    m_spont{index,2} = 0;
                end
            end
%         
%              figure(6)
%              plot(m_spont{index,1})
%              keyboard;
        end
       
    end

    
 %% evoked %%   
    if b_evok
        m_evok = cell(s_numUps,2);
        for index=1:s_numUps
            v_stim = find(v_time>v_evok(index));
            s_stim = v_stim(1);
             if b_total
                s_stim = s_stim - floor(s_baseline/(st_info(1).header.sampleinterval/1000));
            end
            st_infoEvok =st_info;
            st_infoEvok(1).data = v_Vm(s_stim: min(s_stim+s_maxLen,end));
            st_infoEvok(2).data = v_EEG(s_stim: min(s_stim+s_maxLen,end));
            [m_upstates_evok,v_number] = upstate_analysis(st_infoEvok);
            if ~b_total
                m_evok{index,1} = m_upstates_evok{1,1};
                m_evok{index,2} = m_upstates_evok{1,2};
                %0.1 sec before stim and 1 sec after stim
            else  m_evok{index,1} = st_infoEvok(1).data(1:...
                    min(floor((s_baseline + s_trialLen)/(st_info(1).header.sampleinterval/1000)),end));
                  m_evok{index,2} = st_infoEvok(2).data(1:...
                      min(floor((s_baseline + s_trialLen)/(st_info(1).header.sampleinterval/1000)),end));
            end
            
            if ~b_total
                if v_number(1)>length(st_infoEvok(1).data)/10
                    m_evok{index,1} = 0;
                    m_evok{index,2} = 0;
                end
            end     
        end
    end

end