%computes mean time-frequency map for one cell
%mean is computed by cutting all maps to the same length

%ask for number of periods
prompt = {'Number of periods'};
dlg_title = 'Information required';
num_lines = 1;
def = {'1'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
s_nbPeriods = str2num(answer{1});

TFreso = 200;
v_TimeAxisTot = [];
m_TFAbsvm = []; %mean of all intra time-frequency maps
m_TFAbsEEG = [];%mean of all EEG time-frequency maps
s_count = 0;
v_numels=zeros(1,5000); %same number of elements for Vm and EEG
for s_period = 1:s_nbPeriods
    
    str_file = strcat('period', int2str(s_period),'.smr');
    st_info = openSpike2(str_file);
    [m_upstates,v_number,st_info] = upstate_analysis(st_info); 
    
    for s_which = 1:length(v_number)
        
        v_unitVm = m_upstates{s_which,1};
        v_unitEEG = m_upstates{s_which,2};
        b_valid = validity_test(v_unitVm,st_info(1).header.sampleinterval);
        
        if b_valid
%% Compute time-frequency map %%            
            s_count = s_count + 1;
            %mirror before wavelet transform to avoid edge effect
            v_auxVmG = 2*v_unitVm(1)-v_unitVm(end:-1:1);
            v_auxVmR = 2*v_unitVm(end)-v_unitVm(end:-1:1);
            v_unitVm3 = [v_auxVmG,v_unitVm,v_auxVmR];
            v_auxEEGG = 2*v_unitEEG(1)-v_unitEEG(end:-1:1);
            v_auxEEGR = 2*v_unitEEG(end)-v_unitEEG(end:-1:1);
            v_unitEEG3 = [v_auxEEGG,v_unitEEG,v_auxEEGR];
            %compute time-frequency map
            m_vm = f_GaborTransformWait(v_unitVm3,st_info(1).header.Sampling,20,100,TFreso);
            [m_EEG, v_TimeAxis, v_FreqAxis] = ...
            f_GaborTransformWait(v_unitEEG3,st_info(2).header.Sampling,20,100,TFreso);
            %cut mirror parts
            m_vm = m_vm(:,floor(end/3):floor(2*end/3)-1);
            m_EEG = m_EEG(:,floor(end/3):floor(2*end/3)-1);
            v_TimeAxis = v_TimeAxis(1:length(m_vm));
%% compute average %%
            %Intra
            if isempty(m_TFAbsvm)
                m_TFAbsvm = m_vm;
                v_numels(1:size(m_vm,2))=1;
            else s_W = max(size(m_vm,2),size(m_TFAbsvm,2));
                 m_TFAbsvm = fillzeros(m_TFAbsvm,TFreso,s_W) + ...
                            fillzeros(m_vm,TFreso,s_W);
                v_numels(1:size(m_vm,2))= v_numels(1:size(m_vm,2)) + 1;
            end
            %EEG
            if isempty(m_TFAbsEEG)
                m_TFAbsEEG = m_EEG;
            else s_W = max(size(m_EEG,2),size(m_TFAbsEEG,2));
                m_TFAbsEEG = fillzeros(m_TFAbsEEG,TFreso,s_W) + ...
                                fillzeros(m_EEG,TFreso,s_W);
            end
         
            %Time Axis
            if length(v_TimeAxis)> length(v_TimeAxisTot)
                v_TimeAxisTot = v_TimeAxis;
            end
%              
%             figure(12)
%             subplot(1,3,1);
%             f_ImageArray(m_vm,v_TimeAxis, v_FreqAxis,'colormap','jet');
%             subplot(1,3,2);
%             f_ImageArray(m_EEG,v_TimeAxis, v_FreqAxis,'colormap','jet');
%             subplot(2,3,3);
%             plot(v_unitVm)
%             subplot(2,3,6);
%             plot(v_unitEEG);
%             keyboard;
        end
    end
end

   
%Average
for s_bin=1:size(m_TFAbsvm,2)
    m_TFAbsvm(:,s_bin) = m_TFAbsvm(:,s_bin)/v_numels(s_bin);
    m_TFAbsEEG(:,s_bin) = m_TFAbsEEG(:,s_bin)/v_numels(s_bin);
end


figure(4)
Ax(1) = subplot(2,2,1)
    f_ImageArray(m_TFAbsvm, v_TimeAxisTot, v_FreqAxis,'colormap','jet');
    xlabel('Temps (s)');
    ylabel('Fr�quence (Hz)');
    title('Intra')
Ax(2) = subplot(2,2,2)
    f_ImageArray(m_TFAbsEEG, v_TimeAxisTot, v_FreqAxis,'colormap','jet');
    xlabel('Temps (s)');
    ylabel('Fr�quence (Hz)');
    title('EEG')
Ax(3) = subplot(2,2,3)
    plot(0:0.36/1000:0.36/1000*(size(m_TFAbsEEG,2)-1),v_numels(1:size(m_TFAbsEEG,2)));
    title('Nombre de contributions');
Ax(4) = subplot(2,2,4)
    plot(0:0.36/1000:0.36/1000*(size(m_TFAbsEEG,2)-1),v_numels(1:size(m_TFAbsEEG,2)));
    title('Nombre de contributions');
linkaxes(Ax,'x');