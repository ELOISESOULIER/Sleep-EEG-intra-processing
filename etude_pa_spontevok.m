%computes different informations on action potentials for spontaneous and
%evoked upstates

v_filesEvok = [{'c976.4-fqcy.smr'};{'c1036.2_fqcy.smr'};{'c1046.5-fqcy.smr'};...
     {'c1054.2-fqcy.smr'};{'c1307-fqcy.smr'};{'c1520.9-fqcy.smr'};...
     {'c2051_fqcy.smr'};{'c1176.3-evokfqcy1.smr'};{'c1277.6-fqcy2evok.smr'};...
     {'c1615.5-fqcyevok.smr'};{'c2098.4-fqcyEvok.smr'}];
 v_filesSpont = [{'c976.4-fqcy.smr'};{'c1036.2_fqcy.smr'};{'c1046.5-fqcy.smr'};...
     {'c1054.2-fqcy.smr'};{'c1307-fqcy.smr'};{'c1520.9-fqcy.smr'};...
     {'c2051_fqcy.smr'};{'c1176.3-spontfqcy.smr'};{'c1277.6-fqcySpont.smr'};...
     {'c1615.5-fqcyspont.smr'};{'c2098.4-fqcySpont.smr'}];
 
 assert(length(v_filesEvok) == length(v_filesSpont));
 s_nbFiles = length(v_filesEvok);

v_interpaSpont = []; % intervals between action potentials
v_posSpont = []; %position on the upstate
v_seriesSpont = [];  %rank in the upstate
v_valSpont = []; %membrane potential when the action potential is triggered
v_phiSpont = []; %angle of the gamma oscillation of the EEG when the action 
                %potential is triggered

v_interpaEvok = [];  % intervals between action potentials
v_posEvok = [];  %position on the upstate
v_seriesEvok = [];  %rank in the upstate
v_valEvok = []; %membrane potential when the action potential is triggered
v_phiEvok = [];%angle of the gamma oscillation of the EEG when the action 
                %potential is triggered
 
for file = 1:s_nbFiles
    
%% read info and extract upstates    
    info = openSpike2(v_filesEvok{file});
    [~,m_evok,~] = spont_evok_extract(info,false,true);
    info = openSpike2(v_filesSpont{file});
    [m_spont,~,info] = spont_evok_extract(info,true,false);
    Si = info(1).header.sampleinterval;
    
%% Spontaneous %%

    s_len = length(m_spont);
    
    for s_which = 1:s_len
        
        v_unitVmSpont = m_spont{s_which,1};
        v_unitEEGSpont = m_spont{s_which,2};
        
        v_diff = (v_unitVmSpont(2:end) - v_unitVmSpont(1:end-1))/360*1000;
        %test length and presence of action potentials
        b_len = length(v_unitVmSpont)*Si/1000 > 300;
        b_pa = sum(v_diff >10) >1.5; %at least 2 action potentials on the upstate
        
        if b_pa && b_len  
            
            %find the action potentials
            v_Index = find(v_diff >10)'; 
            v_swift =[0;v_Index];
            v_aux = ((v_swift(2:end) - v_swift(1:end-1))>2)'; %not twice the same
            v_Index = v_Index(v_aux);
            
            %store position, value, interval with previous one
            v_posSpont = [v_posSpont; v_Index./length(v_unitVmSpont)];
            v_valSpont = [v_valSpont;v_unitVmSpont(v_Index)'];
            v_interpaSpont = [v_interpaSpont;(v_Index(2:end) - v_Index(1:end-1))];
            v_seriesSpont = [v_seriesSpont;s_which*ones(length(v_Index)-1,1)];
            
            % instantaneous value of phase of the gamma oscillation
            v_FiltGamma	= f_DesignIIRfilter(info(2).header.Sampling,[30 100],[29 101],[1,100]);
            v_filtGamEEG = f_FilterIIR(v_unitEEGSpont,v_FiltGamma);
            v_analyticGamma = hilbert(v_filtGamEEG);
            s_phase = angle(v_analyticGamma);
            v_phiSpont = [v_phiSpont;s_phase(v_Index)];
            
            %visualize phase of gamma vs signal
            figure(18)
            Ax(1)= subplot(2,1,1)
                plot(v_unitVmSpont);
                title('Signal brut')
                xlabel('Temps (pts)')
                ylabel('Potentiel (mV)')
            Ax(2) = subplot(2,1,2);
                plot(s_phase);
                hold on;
                plot(v_Index,s_phase(v_Index),'*');
                title('Phase correspondante des gamma de l EEG');
                xlabel('Temps (pts)');
                ylabel('Phase (rad)');
                hold off;
            linkaxes(Ax,'x');
            keyboard;
        end
    end

    v_interpaSpont = v_interpaSpont*Si*10^-3;

%% Evoked %%

    s_len = length(m_evok);
    
    for s_which=1:s_len
        
        v_unitVmEvok= m_evok{s_which,1};
        v_unitEEGEvok = m_evok{s_which,2};
        
        v_diff = (v_unitVmEvok(2:end) - v_unitVmEvok(1:end-1))/360*1000;
        %test length and presence of action potentials
        b_len = length(v_unitVmEvok)*Si/1000 > 300;
        b_pa = sum(v_diff >10) >1.5;%at least 2 action potentials on the upstate
        
        if b_pa && b_len  
            
            %find action potentials
            v_Index = find(v_diff >10)'; 
            v_swift =[0;v_Index];
            v_aux = ((v_swift(2:end) - v_swift(1:end-1))>2)'; %not twice the same
            v_Index = v_Index(v_aux);
            
            %store position, value, interval with previous action potential
            v_posEvok = [v_posEvok; v_Index./length(v_unitVmEvok)];
            v_valEvok = [v_valEvok;v_unitVmEvok(v_Index)'];
            v_interpaEvok = [v_interpaEvok;(v_Index(2:end) - v_Index(1:end-1))];
            v_seriesEvok = [v_seriesEvok;s_which*ones(length(v_Index)-1,1)];
            
            %instantaneous phase of gamma oscillation
            v_FiltGamma	= f_DesignIIRfilter(info(2).header.Sampling,[30 100],[29 101],[1,100]);
            v_filtGamEEG = f_FilterIIR(v_unitEEGEvok,v_FiltGamma);
            v_analyticGamma = hilbert(v_filtGamEEG);
            s_phase = angle(v_analyticGamma);
            v_phiEvok= [v_phiEvok;s_phase(v_Index)];
            
            % visualize phase of gamma vs signal
            figure(18)
            Ax(1)=subplot(2,1,1)
                plot(v_unitVmEvok);
                title('Signal brut')
                xlabel('Temps (pts)')
                ylabel('Potentiel (mV)')
            Ax(2)=subplot(2,1,2);
                plot(s_phase);
                hold on;
                plot(v_Index,s_phase(v_Index),'*');
                title('Phase correspondante des gamma de l EEG');
                xlabel('Temps (pts)');
                ylabel('Phase (rad)');
                hold off;
            linkaxes(Ax,'x');
        end
    end

    v_interpaEvok = v_interpaEvok*Si*10^-3;

%% Visualize %%
    
    figure(12)
        histogram(v_interpaSpont,40);
        hold on;
        histogram(v_interpaEvok,40);
        hold off;
        legend('Spontan�','Evoqu�')
        xlabel('time, ms')
        ylabel('nb of intervals')
        title('Distribution des intervalles inter pa')
    figure(13);
        histogram(v_posSpont,40);
        hold on;
        histogram(v_posEvok,40);
        hold off;
        legend('Spontan�','Evoqu�')
        xlabel('Proportion de lupstate')
        title('Distribution des positions des pa')
    figure(14)
        histogram(v_phiSpont,40);
        hold on;
        histogram(v_phiEvok,40);
        hold off;
        legend('Spontan�','Evoqu�')
        xlabel('Valeur de la phase')
        title('Position des pa par rapport � la phase des oscillations gamma de lEEG')

    %figure(14)
    %plot(v_position,v_inter_pa,'*');
    %xlabel('latence, ms')
    %ylabel('intervalle inter pa, ms')
    %title('intervalle avec le prochain pa en fonction de la latence')

    %n en fonction de n-1 : prendre en compte les trous !
    %x=[];%n-1
    %y=[]; %n
    %for i=1:length(v_series)-1
    %    for k = 3:v_series(i)
    %        x = [x,v_inter_pa(sum(v_series(1:i)) +k-2)];
    %        y = [y,v_inter_pa(sum(v_series(1:i)) +k-1)];
    %    end
    %end
    %figure()
    %plot(x,y,'*');

end
