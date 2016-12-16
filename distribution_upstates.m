%% computes distribution of values on the upstates of all recordings of one cell

clear
%distribution des valeurs des upstates
%pour r�pondre � deux questions :
%- la question de St�phane de savoir quelle est la distribution des valeurs
%apr�s d�coupage des upstates par rapport � la distribution initiale
%- la question de savoir si la ditribution des seuils de potentiel de
%l'intra � l'apparition d'un potentiel d'action est comparable � la
%distribution des valeurs de l'upstate
% 
nb_periods = 7;
v_values=[];
v_upValues=[];
s_count=0;

for s_period=1:nb_periods
    
    file = strcat('period', int2str(s_period),'.smr');
    info=openSpike2(file);
    [m_upstates,number,info] =upstate_analysis(info);
    v_values = [v_values,info(1).data];
    for index=1:length(number)
        s_count=s_count+1;
        v_upValues = [v_upValues,m_upstates{index,1}];
    end
    
end

figure(40)
histogram(v_values,100)
hold on
histogram(v_upValues,100)
hold off
title('Distribution des valeurs des upstates apr�s d�coupage')
xlabel('Potentiel (mV)')

    
% figure()
% [v_N,v_edges]=histcounts(v_upValues);
% v_N = v_N/length(v_upValues);
% bar(v_edges(1:end-1),v_N);
% hold on;
% [v_N,v_edges]=histcounts(Vm);
% v_N = v_N/2/length(Vm);
% bar(v_edges(1:end-1),v_N,'r');
% hold off;
% legend('Distribution totale','Au declenchement du pa')
% title('Potentiel au declenchement d un potentiel d action')
% xlabel('Potentiel (mV)')

% info=openSpike2('c2051_fqcy.smr')
% [~,m_evok, ~]=spont_evok_extract(info,false,true);
% info=openSpike2('c2051_fqcy.smr')
% [m_spont,~,info]=spont_evok_extract(info,true,false);
% 
% v_upValuesSpont = [];
% v_upValuesEvok = [];
% v_EEGValuesSpont=[];
% v_EEGValuesEvok =[];
% s_counts=0;
% s_counte=0;
% 
% for index=1:length(m_spont)
%     s_counts=s_counts+1;
%     v_upValuesSpont = [v_upValuesSpont,m_spont{index,1}];
%     v_EEGValuesSpont = [v_EEGValuesSpont,m_spont{index,2}];
%     figure(18)
%     time = 0:info(1).header.sampleinterval/1000:info(1).header.sampleinterval/1000*(length(m_spont{index,1})-1);
%     subplot(2,2,1)
%         plot(time,m_spont{index,1},'b')
%         title('Intra');xlabel('Temps (ms)'); ylabel('Potentiel (mV)');
%         hold on;
%      subplot(2,2,2)
%         plot(time,m_spont{index,2},'b')
%         title('EEG');xlabel('Temps (ms)'); ylabel('Potentiel (uV)');
%         hold on;
%     
% end
% for index = 1:length(m_evok)
%     s_counte = s_counte+1;
%     v_upValuesEvok = [v_upValuesEvok,m_evok{index,1}];
%     v_EEGValuesEvok = [v_EEGValuesEvok, m_evok{index,2}];
%     time = 0:info(1).header.sampleinterval/1000:info(1).header.sampleinterval/1000*(length(m_evok{index,1})-1);
%     figure(18)
%     subplot(2,2,1)
%         plot(time,m_evok{index,1},'r')
%         hold on;
%     subplot(2,2,2)
%         plot(time,m_evok{index,2},'r')
%         hold on;
%     
% end
% 
% figure(18)
% subplot(2,2,3)
%     histogram(v_upValuesSpont)
%     hold on;
%     histogram(v_upValuesEvok)
%     hold off;
%     title('Distribution potentiel des upstates');xlabel('Potentiel (mV)')
%     legend('Spontan�','Evoqu�');
% subplot(2,2,4)
%     histogram(v_EEGValuesSpont)
%     hold on;
%     histogram(v_EEGValuesEvok)
%     hold off;
%     title('Distribution potentiel des upstates');xlabel('Potentiel(uV)')
%     legend('Spontan�','Evoqu�');
%     
%     
% values =zeros(12,1);
% values(1) = mean(v_upValuesEvok);
% values(2) = median(v_upValuesEvok);
% values(3) = std(v_upValuesEvok);
% values(4) = mean(v_upValuesSpont);
% values(5) = median(v_upValuesSpont);
% values(6) = std(v_upValuesSpont);
% 
% values(7) = mean(v_EEGValuesEvok);
% values(8) = median(v_EEGValuesEvok);
% values(9) = std(v_EEGValuesEvok);
% values(10) = mean(v_EEGValuesSpont);
% values(11) = median(v_EEGValuesSpont);
% values(12) = std(v_EEGValuesSpont);
