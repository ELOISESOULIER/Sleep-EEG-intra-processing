% script used to do the analysis of a particular trace

addpath('/home/eloise/Stages/Rythm-ICM/Spike2-data/')
addpath(genpath('/home/eloise/Stages/Rythm-ICM/Matlab/Toolbox/'))

params = getParams();

info = openSpike2('c1520.smr');
v_Vm = info(1).data*100;
v_EEG = info(2).data;

[m_TFbefore,v_TimeAxis,v_FreqAxis] = compute_TF(v_Vm,info(1).header.Sampling,...
             10,100,params.s_TFreso);
figure(8)
plot(v_TimeAxis,v_Vm);
hold on;
         
v_Vm = traitementpa_general(v_Vm,30,0.36);

figure(8)
plot(v_TimeAxis,v_Vm)
xlabel('Temps (s)');ylabel('Potentiel (mV)'); title('Pa supprimés');
hold off;

m_TFafter = compute_TF(v_Vm,info(1).header.Sampling,...
             9,100,params.s_TFreso);

m_TFEEG = compute_TF(v_EEG,info(2).header.Sampling,9,100,params.s_TFreso);    
     
 figure(50)
 subplot(2,1,1)
 v_clim = [prctile(m_TFafter(:),1),prctile(m_TFafter(:),96)];
 f_ImageArray(m_TFafter,v_TimeAxis,v_FreqAxis,'colormap','jet','limits',v_clim);
 title('Vm'); ylabel('Hz'); xlabel('s');
 colorbar
 subplot(2,1,2)
 v_clim = [prctile(m_TFEEG(:),1),prctile(m_TFEEG(:),98)];
 f_ImageArray(m_TFEEG,v_TimeAxis(1:end-1),v_FreqAxis,'colormap','jet','limits',v_clim);
 title('EEG');ylabel('Hz'); xlabel('s');
 colorbar
 
 v_partialVm1 = v_Vm(floor(1690/0.36):floor(2350/0.36));
 v_partialVm2 = v_Vm(floor(4440/0.36):floor(5390/0.36));
 
 v_partialEEG1 = v_EEG(floor(1690/0.36):floor(2350/0.36));
 v_partialEEG2 = v_EEG(floor(4440/0.36):floor(5390/0.36));
 
 m_partialVm1 = m_TFafter(:,floor(1690/0.36):floor(2350/0.36));
 m_partialVm2 = m_TFafter(:,floor(4440/0.36):floor(5390/0.36));
 
 m_partialEEG1 = m_TFEEG(:,floor(1690/0.36):floor(2350/0.36));
 m_partialEEG2 = m_TFEEG(:,floor(4440/0.36):floor(5390/0.36));
 
 v_partialTime1 = v_TimeAxis(floor(1690/0.36):floor(2350/0.36));
 v_partialTime2 = v_TimeAxis(floor(4440/0.36):floor(5390/0.36));
 
 %% Visualize
 
 figure(1)
 plot(v_partialTime1,v_partialVm1);
 xlabel('Temps (s)');title('Vm 1');
 
 figure(2)
 v_clim = [prctile(m_partialVm1(:),0),prctile(m_partialVm1(:),99.5)];
 f_ImageArray(m_partialVm1,v_partialTime1,v_FreqAxis,'colormap','jet','limits',v_clim);
 xlabel('Temps (s)');title('Vm 1');
 colorbar
 
 figure(3)
 plot(v_partialTime2,v_partialVm2);
 xlabel('Temps (s)');title('Vm 2');
 
 figure(4)
 v_clim = [prctile(m_partialVm2(:),0),prctile(m_partialVm2(:),99.5)];
 f_ImageArray(m_partialVm2,v_partialTime2,v_FreqAxis,'colormap','jet','limits',v_clim);
 xlabel('Temps (s)');title('Vm 2');
 colorbar
 
 figure(5)
 plot(v_partialTime1,v_partialEEG1);
 xlabel('Temps (s)');title('EEG 1');
 
 figure(6)
 v_clim = [prctile(m_partialEEG1(:),1),prctile(m_partialEEG1(:),99)];
 f_ImageArray(m_partialEEG1,v_partialTime1,v_FreqAxis,'colormap','jet','limits',v_clim);
 xlabel('Temps (s)');title('EEG 1');
 colorbar
 
 figure(7) 
 plot(v_partialTime2,v_partialEEG2);
 xlabel('Temps (s)');title('EEG 2');
 
 figure(8)
 v_clim = [prctile(m_partialEEG2(:),1),prctile(m_partialEEG2(:),99)];
 f_ImageArray(m_partialEEG2,v_partialTime2,v_FreqAxis,'colormap','jet','limits',v_clim);
 xlabel('Temps (s)');title('EEG 2');
 colorbar