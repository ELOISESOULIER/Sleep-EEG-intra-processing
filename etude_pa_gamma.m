%computes different informations on action potentials

s_nbPeriods = 7;
v_inter_pa = []; %intervals between action potentials
v_position = []; %position on the upstate
v_series = []; %rank in the upstate
v_values = []; %membrane potential when the action potential is triggered
v_phases = []; %angle of the gamma oscillation of the EEG when the action 
                %potential is triggered
    
for s_period=1:s_nbPeriods
    
    fprintf('period %1f\n',s_period);
    
%% read file and extract upstates
    str_file = strcat('period', int2str(s_period),'.smr');
    st_info = openSpike2(str_file);
    [m_upstates,~,st_info] = upstate_analysis(st_info);
    s_SI = st_info(1).header.sampleinterval;

    for s_which = 1:length(m_upstates)
        
        v_Vm = m_upstates{s_which,1};
        v_EEG = m_upstates{s_which,2};
        
        v_diffVm = (v_Vm(2:end) - v_Vm(1:end-1))/360*1000;
        %test length and presence of action potentials
        b_len = length(v_Vm)*s_SI/1000 > 300; 
        b_pa = sum(v_diffVm >10) >1.5; %at least 2
        
        if b_pa && b_len   
            
            %find action potentials
            v_Index = find(v_diffVm >10)'; 
            v_swift = [0;v_Index];
            v_aux = ((v_swift(2:end) - v_swift(1:end-1))>2)';
            
            %store position, value, interval with previous action potential
            v_Index = v_Index(v_aux);
            v_position = [v_position; v_Index./length(v_Vm)];
            v_values = [v_values;v_Vm(v_Index)'];
            v_inter_pa = [v_inter_pa;(v_Index(2:end) - v_Index(1:end-1))];
            v_series = [v_series;s_which*ones(length(v_Index)-1,1)];
            %instantaneous phase of gamma oscillation
            v_FiltGam	= f_DesignIIRfilter(st_info(2).header.Sampling,...
                [30 100],[29 101],[1,100]);
            v_filtGamEEG = f_FilterIIR(v_EEG,v_FiltGam);
            v_analyticGamma = hilbert(v_filtGamEEG);
            s_phase = angle(v_analyticGamma);
            v_phases = [v_phases;s_phase(v_Index)];
            
            %visualize phase of gamma vs signal
            figure(18)
            time = 0: s_SI/1000:s_SI/1000*(length(v_Vm)-1);
            Ax(1)=subplot(2,1,1)
                plot(time,v_Vm);
                title('Signal brut')
                xlabel('Temps (ms)')
                ylabel('Potentiel (mV)')
            Ax(2)=subplot(2,1,2);
                plot(time,s_phase');
                hold on;
                plot(v_Index*s_SI/1000,s_phase(v_Index),'*');
                title('Phase correspondante des gamma de l EEG');
                xlabel('Temps (ms)');
                ylabel('Phase (rad)');
            hold off;
            linkaxes(Ax,'x');
            keyboard;
          
       
        end

    end
end
v_inter_pa = v_inter_pa*s_SI*10^-3; %convert from points to duration

%% visualize %%

figure(12)
histogram(v_inter_pa,50)
xlabel('time, ms')
ylabel('nb of intervals')
title('Distribution des intervalles inter pa')

figure(13);
histogram(v_position,50);
xlabel('Proportion de lupstate')
title('Distribution des positions des pa')

figure(14)
histogram(v_phases,50);
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
        
        
        