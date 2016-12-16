%computes mean power spectrum on one period for gamma and spindle
%frequencies
%input can be a string which is interpreted as the name of the spike2 file,
%or a cell, which is considered as the output of the upstate_analysis
%function

function [v_meanSpVm,v_meanGamVm,v_meanSpEEG,v_meanGamEEG] = mean_spec(str_file)


if ischar(str_file)
    
    st_info = openSpike2(str_file);
    [m_upstates,~,st_info] = upstate_analysis(st_info);
    s_SI = st_info(1).header.sampleinterval;
    s_sampling = st_info(1).header.Sampling;
    
else
    m_upstates = str_file;
    s_SI = 360;
    s_sampling = 2.7778e+03;
    
end

%% compute power spectrum for each frequency interval %%
m_powSpectra = repmat(struct('sp_vm',[],'gam_vm',[],...
    'sp_EEG',[],'gam_EEG',[]),length(m_upstates),1);

s_valid = 0;
for s_up = 1:length(m_upstates)
    
     % check length and absence of action potentials
     [bool,v_newVm] =validity_test(m_upstates{s_up,1},s_SI);
     
     if bool
        s_valid = s_valid + 1;
         [m_powSpectra(s_up).sp_vm, m_powSpectra(s_up).gam_vm,...
            m_powSpectra(s_up).sp_EEG, m_powSpectra(s_up).gam_EEG]...
             = filtre_spectre(v_newVm,m_upstates{s_up,2},s_sampling);
         
     end
    
end

fprintf('number of valid upstates for this period %f\n',s_valid);

%%  compute mean spectrum over one period for both frequency bands

v_meanSpVm = [];
for s_up=1:length(m_upstates)
    if ~isempty(m_powSpectra(s_up).sp_vm)
        if isempty(v_meanSpVm)
            v_meanSpVm = m_powSpectra(s_up).sp_vm;
        else
            N = min(length(v_meanSpVm),length(m_powSpectra(s_up).sp_vm));
            v_meanSpVm = v_meanSpVm(1:N);
            m_powSpectra(s_up).sp_vm = m_powSpectra(s_up).sp_vm(1:N);
            v_meanSpVm = (v_meanSpVm+m_powSpectra(s_up).sp_vm)/2;
        end
    end
end

v_meanGamVm = [];
for s_up=1:length(m_upstates)
    if ~isempty(m_powSpectra(s_up).gam_vm)
        if isempty(v_meanGamVm)
            v_meanGamVm = m_powSpectra(s_up).gam_vm;
        else
            N = min(length(v_meanGamVm),length(m_powSpectra(s_up).gam_vm));
            v_meanGamVm = v_meanGamVm(1:N);
            m_powSpectra(s_up).gam_vm = m_powSpectra(s_up).gam_vm(1:N);
            v_meanGamVm = (v_meanGamVm+m_powSpectra(s_up).gam_vm)/2;
        end
    end
end

v_meanSpEEG = [];
for s_up=1:length(m_upstates)
    if ~isempty(m_powSpectra(s_up).sp_EEG)
        if isempty(v_meanSpEEG)
            v_meanSpEEG = m_powSpectra(s_up).sp_EEG;
        else
            N = min(length(v_meanSpEEG),length(m_powSpectra(s_up).sp_EEG));
            v_meanSpEEG = v_meanSpEEG(1:N);
            m_powSpectra(s_up).sp_EEG = m_powSpectra(s_up).sp_EEG(1:N);
            v_meanSpEEG = (v_meanSpEEG+m_powSpectra(s_up).sp_EEG)/2;
        end
    end
end

v_meanGamEEG = [];
for s_up=1:length(m_upstates)
    if ~isempty(m_powSpectra(s_up).gam_EEG)
        if isempty(v_meanGamEEG)
            v_meanGamEEG = m_powSpectra(s_up).gam_EEG;
        else
            N = min(length(v_meanGamEEG),length(m_powSpectra(s_up).gam_EEG));
            v_meanGamEEG = v_meanGamEEG(1:N);
            m_powSpectra(s_up).gam_EEG = m_powSpectra(s_up).gam_EEG(1:N);
            v_meanGamEEG = (v_meanGamEEG+m_powSpectra(s_up).gam_EEG)/2;
        end
    end
end

end
