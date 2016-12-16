s_nbPeriod = 3;
s_freqMin = 20;
s_freqMax = 60;
s_corrBin = 4; % correlation interval width
s_nbBands = 2*(s_freqMax-s_freqMin);
v = s_freqMin+1:s_corrBin/2:s_freqMax;

m_corrstot =[];
s_count = 0;
for s_period = 1:s_nbPeriod
    file = strcat('period', int2str(s_period),'.smr');
    st_info = openSpike2(file);
    [m_upstates,number,st_info] = upstate_analysis(st_info); 
%% select "good" upstates and make time-frequency map %%    
    for s_which = 1:length(number)
        v_unitVm = m_upstates{s_which,1};
        v_unitEEG = m_upstates{s_which,2};
        %test length and absence of action potential
        b_valid = validity_test(v_unitVm,st_info(1).header.sampleinterval);
        if b_valid
            s_count = s_count + 1;
            [m_vm, v_TimeAxis, v_FreqAxis] = ...
            compute_TF(v_unitVm,st_info(1).header.Sampling,...
                                    s_freqMin,s_freqMax,s_nbBands);
            m_EEG = compute_TF(v_unitEEG,st_info(2).header.Sampling...
                                            ,s_freqMin,s_freqMax,s_nbBands);
%% Build inter-frequency interval correlation matrix %%                                        
            m_corrsAux = zeros(s_nbBands); %one2one correlations
            for s_abs = 1:s_nbBands
                v_vmref = m_vm(s_abs,:) - mean(m_vm(s_abs,:));
                for s_ord = 1:s_nbBands
                    m_corrsAux(s_abs,s_ord) = ...
                     max(xcorr(v_vmref,m_EEG(s_ord,:)-mean(m_EEG(s_ord,:)),'coeff'));
                end
            end
            m_corrs = zeros(floor(s_nbBands/s_corrBin));%reduction
            for s_abs = 1:floor(s_nbBands/s_corrBin)
                for s_ord = 1:floor(s_nbBands/s_corrBin)
                    coeff = [];
                    for index = 1:s_corrBin
                        coeff = [coeff,m_corrsAux((s_abs-1)*s_corrBin+index,(s_ord-1)*s_corrBin+index)]
                    end
                    m_corrs(s_abs,s_ord) = mean(coeff);
                end
            end
                    
            if isempty(m_corrstot) %mean of all correlation matrices
                m_corrstot = m_corrs;
            else m_corrstot = m_corrstot+m_corrs;
            end
       
%% visualize %%
            time = 0:st_info(1).header.sampleinterval/1000:st_info(1).header.sampleinterval/1000*(length(v_unitVm)-1);
            figure(7)
                subplot(3,2,1)
                plot(time,v_unitVm);
                title('Intra');xlabel('Temps (sec)');ylabel('Fr�quence (Hz)');
            subplot(3,2,2)
                plot(time,v_unitEEG)
                title('EEG');xlabel('Temps (sec)'); ylabel('Fr�quence (Hz)');
            subplot(3,2,3)
                f_ImageArray(m_vm,v_TimeAxis, v_FreqAxis,'colormap','jet');
                title('Intra');
                xlabel('Temps (sec)'); ylabel('Fr�quence (Hz)');
            subplot(3,2,4);
                f_ImageArray(m_EEG,v_TimeAxis, v_FreqAxis,'colormap','jet')
                title('EEG');
                xlabel('Temps (sec)'); ylabel('Fr�quence (Hz)');
            subplot(3,1,3);
                f_ImageArray(m_corrs,v,v,'colormap','jet');
                title('Correlations inter bandes de fr�quence');
                xlabel('Fr�quence Intra (Hz)'); ylabel('Fr�quence EEG (Hz)');            
        
          
            figure(17)
            f_ImageArray(m_corrstot,v,v,'colormap','jet');
            title('Correlations inter bandes de fr�quence');
            xlabel('Fr�quence Intra (Hz)'); ylabel('Fr�quence EEG (Hz)');  

        end
        
     end
end
                    
              
m_corrstot=m_corrstot./s_count;

figure(2)
f_ImageArray(m_corrstot,v,v,'colormap','jet');
title('Correlation moyenne entre bandes de fr�quence');
xlabel('Fr�quence Intra (Hz)'); ylabel('Fr�quence EEG (Hz)');  