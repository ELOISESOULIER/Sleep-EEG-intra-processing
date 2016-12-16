% Extract power from mean TF maps saved as .mat

v_cells = [{'976'};{'1036'};{'1045'};{'1054'};{'1277'};...
    {'1307'};{'1520'};{'1615'};{'2051'};{'2098'}]; %{'1176'};

params = getParams();
s_TFreso = params.s_TFreso;
s_sampInt = 360;
v_FreqAxis = linspace(100,8,s_TFreso);
s_cut = floor(150/s_sampInt*1000);

m_TFtotSpont = [];
m_TFtotEvok = [];
m_TFtotSpontNorm = [];
m_TFtotEvokNorm = [];

v_meanSpLoSpont = [];
v_meanSpHiSpont = [];
v_meanGamLoSpont = [];
v_meanGamMidSpont = [];
v_meanGamHiSpont = [];

v_meanSpLoEvok = [];
v_meanSpHiEvok = [];
v_meanGamLoEvok = [];
v_meanGamMidEvok = [];
v_meanGamHiEvok = [];

for s_cellNum = 1:length(v_cells)

    load(['/home/eloise/Stages/Rythm-ICM/Matlab/' v_cells{s_cellNum} '-trials/totalSpontNorm']);
    load(['/home/eloise/Stages/Rythm-ICM/Matlab/' v_cells{s_cellNum} '-trials/totalEvokNorm']);
    
    load(['/home/eloise/Stages/Rythm-ICM/Matlab/' v_cells{s_cellNum} '-trials/totalSpont']);
    load(['/home/eloise/Stages/Rythm-ICM/Matlab/' v_cells{s_cellNum} '-trials/totalEvok']);
    
    if isempty(m_TFtotSpont)
        m_TFtotSpont = m_TFAbsEEGSpont;
    else
        m_TFtotSpont = m_TFtotSpont + m_TFAbsEEGSpont(:,1:length(m_TFtotSpont));
    end
    
    if isempty(m_TFtotSpontNorm)
        m_TFtotSpontNorm = m_TFAbsEEGSpontNorm;
    else
        m_TFtotSpontNorm = m_TFtotSpontNorm + m_TFAbsEEGSpontNorm(:,1:length(m_TFtotSpontNorm));
    end
    
    if isempty(m_TFtotEvok)
        m_TFtotEvok = m_TFAbsEEGEvok;
    else
        m_TFtotEvok = m_TFtotEvok + m_TFAbsEEGEvok(:,1:length(m_TFtotEvok));
    end
    
    if isempty(m_TFtotEvokNorm)
        m_TFtotEvokNorm = m_TFAbsEEGEvokNorm;
    else
        m_TFtotEvokNorm = m_TFtotEvokNorm + m_TFAbsEEGEvokNorm(:,1:length(m_TFtotEvokNorm));
    end
    
    v_TimeAxis = 0:s_sampInt/1000000:s_sampInt/1000000*(length(m_TFtotSpont)-1);

    [v_spLoSpont,v_spLoEvok,v_spHiSpont,v_spHiEvok,v_gamLoSpont,v_gamLoEvok,...
    v_gamMidSpont,v_gamMidEvok,v_gamHiSpont,v_gamHiEvok] = ...
    power_extract(m_TFAbsEEGSpontNorm(:,s_cut:end),...
                    m_TFAbsEEGEvokNorm(:,s_cut:end),v_FreqAxis);
                
    if isempty(v_meanSpLoSpont)
        v_meanSpLoSpont = v_spLoSpont;
        v_meanSpHiSpont = v_spHiSpont;
        v_meanGamLoSpont = v_gamLoSpont;
        v_meanGamMidSpont = v_gamMidSpont;
        v_meanGamHiSpont = v_gamHiSpont;

    else v_meanSpLoSpont = v_meanSpLoSpont + v_spLoSpont(1:length(v_meanSpLoSpont));
         v_meanSpHiSpont = v_meanSpHiSpont + v_spHiSpont(1:length(v_meanSpHiSpont));
         v_meanGamLoSpont = v_meanGamLoSpont + v_gamLoSpont(1:length(v_meanGamLoSpont));
         v_meanGamMidSpont = v_meanGamMidSpont + v_gamMidSpont(1:length(v_meanGamMidSpont));
         v_meanGamHiSpont = v_meanGamHiSpont + v_gamHiSpont(1:length(v_meanGamHiSpont));

    end
    
    if isempty(v_meanSpLoEvok)
        v_meanSpLoEvok = v_spLoEvok;
        v_meanSpHiEvok = v_spHiEvok;
        v_meanGamLoEvok = v_gamLoEvok;
        v_meanGamMidEvok = v_gamMidEvok;
        v_meanGamHiEvok = v_gamHiEvok;

    else
        v_meanSpLoEvok = v_meanSpLoEvok + v_spLoEvok(1:length(v_meanSpLoEvok));
        v_meanSpHiEvok = v_meanSpHiEvok + v_spHiEvok(1:length(v_meanSpHiEvok));
        v_meanGamLoEvok = v_meanGamLoEvok + v_gamLoEvok(1:length(v_meanGamLoEvok));
        v_meanGamMidEvok = v_meanGamMidEvok + v_gamMidEvok(1:length(v_meanGamMidEvok));
        v_meanGamHiEvok = v_meanGamHiEvok + v_gamHiEvok(1:length(v_meanGamHiEvok));

    end
    
    figure(3)
    subplot(1,2,1)
        f_ImageArray(m_TFtotSpont,v_TimeAxis,v_FreqAxis,'colormap','jet');
        title('Spontané');
    subplot(1,2,2)
        f_ImageArray(m_TFtotEvok,v_TimeAxis,v_FreqAxis,'colormap','jet');
        title('Évoqué')

    figure(4)
    subplot(1,2,1)
        v_clim = [prctile(m_TFtotSpontNorm(:),1),prctile(m_TFtotSpontNorm(:),100)];
        f_ImageArray(m_TFtotSpontNorm,v_TimeAxis,v_FreqAxis,'colormap','jet','limits',v_clim);
        title('Spontané');
    subplot(1,2,2)
        v_clim = [prctile(m_TFtotEvokNorm(:),1),prctile(m_TFtotEvokNorm(:),100)];
        f_ImageArray(m_TFtotEvokNorm,v_TimeAxis,v_FreqAxis,'colormap','jet','limits',v_clim);
        title('Évoqué');
      
    figure(34)
    subplot(5,1,1)
        plot(v_TimeAxis,v_spLoSpont,'b');hold on;
        plot(v_TimeAxis,v_spLoEvok,'r');%hold off;
        legend('Spontané','Évoqué');title('Spindle 9-12');
        xlabel('Temps (ms)'); ylabel('Puissance (muV^2)');
    subplot(5,1,2)
        plot(v_TimeAxis,v_spHiSpont,'b');hold on;
        plot(v_TimeAxis,v_spHiEvok,'r');%hold off;
        title('Spindle 12-16');
        xlabel('Temps (ms)'); ylabel('Puissance (uV^2)');
     subplot(5,1,3)
        plot(v_TimeAxis,v_gamLoSpont,'b');hold on;
        plot(v_TimeAxis,v_gamLoEvok,'r');%hold off;
        title('Gamma 30-60');
        xlabel('Temps (ms)'); ylabel('Puissance (uV^2)');
     subplot(5,1,4)
        plot(v_TimeAxis,v_gamMidSpont,'b');hold on;
        plot(v_TimeAxis,v_gamMidEvok,'r');%hold off;
        title('Gamma 60-80');
        xlabel('Temps (ms)'); ylabel('Puissance (uV^2)');  
    subplot(5,1,5)
        plot(v_TimeAxis,v_gamHiSpont,'b');hold on;
        plot(v_TimeAxis,v_gamHiEvok,'r');%hold off;
        title('Gamma 80-100');
        xlabel('Temps (ms)'); ylabel('Puissance (uV^2)');   
        
        saveas(gcf,['/home/eloise/Stages/Rythm-ICM/Matlab/' v_cells{s_cellNum} '-trials/power-extraction.png']);
        
end


figure(3)
subplot(1,2,1)
    f_ImageArray(m_TFtotSpont,v_TimeAxis,v_FreqAxis,'colormap','jet');
    title('Spontané');
subplot(1,2,2)
    f_ImageArray(m_TFtotEvok,v_TimeAxis,v_FreqAxis,'colormap','jet');
    title('Évoqué')

figure(4)
subplot(1,2,1)
    v_clim = [prctile(m_TFtotSpontNorm(:),1),prctile(m_TFtotSpontNorm(:),100)];
    f_ImageArray(m_TFtotSpontNorm,v_TimeAxis,v_FreqAxis,'colormap','jet','limits',v_clim);
    title('Spontané');
subplot(1,2,2)
    v_clim = [prctile(m_TFtotEvokNorm(:),1),prctile(m_TFtotEvokNorm(:),100)];
    f_ImageArray(m_TFtotEvokNorm,v_TimeAxis,v_FreqAxis,'colormap','jet','limits',v_clim);
    title('Évoqué');

  v_TimeAxis = v_TimeAxis(s_cut:end);

figure(34)
subplot(5,1,1)
    plot(v_TimeAxis,v_meanSpLoSpont,'b');hold on;
    plot(v_TimeAxis,v_meanSpLoEvok,'r');%hold off;
    legend('Spontané','Évoqué');title('Spindle 9-12');
    xlabel('Temps (ms)'); ylabel('Puissance (muV^2)');
subplot(5,1,2)
    plot(v_TimeAxis,v_meanSpHiSpont,'b');hold on;
    plot(v_TimeAxis,v_meanSpHiEvok,'r');%hold off;
    title('Spindle 12-16');
    xlabel('Temps (ms)'); ylabel('Puissance (uV^2)');
 subplot(5,1,3)
    plot(v_TimeAxis,v_meanGamLoSpont,'b');hold on;
    plot(v_TimeAxis,v_meanGamLoEvok,'r');%hold off;
    title('Gamma 30-60');
    xlabel('Temps (ms)'); ylabel('Puissance (uV^2)');
 subplot(5,1,4)
    plot(v_TimeAxis,v_meanGamMidSpont,'b');hold on;
    plot(v_TimeAxis,v_meanGamMidEvok,'r');%hold off;
    title('Gamma 60-80');
    xlabel('Temps (ms)'); ylabel('Puissance (uV^2)');  
subplot(5,1,5)
    plot(v_TimeAxis,v_meanGamHiSpont,'b');hold on;
    plot(v_TimeAxis,v_meanGamHiEvok,'r');%hold off;
    title('Gamma 80-100');
    xlabel('Temps (ms)'); ylabel('Puissance (uV^2)');   



