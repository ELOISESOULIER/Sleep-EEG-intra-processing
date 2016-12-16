%Correspondances numéros - graphes :
%1: M1L
%2: S1L
%3: S1R
%4: S1L & S1R
%5: S1L & S1R & M1R

%% Chemins vers les données et vers la librairie SON
s_path = '/home/eloise/Stages/Rythm-ICM/Adrien - All WoD';
Toolbox_path = '/home/eloise/Stages/Rythm-ICM/Matlab/Toolbox/';

%% Paramètres

%numéro du canal qui contient les données
whereisS1L = 3;
whereisM1L = 4;
whereisS1R = 5;

%intervalle d'échantillonage objectif (si tous ne sont pas échantillonés de
%la meme façon)
s_SI = 6600;

%longueur minimale du signal en points (pour ramener tous les signaux à la
%meme longueur)
s_minLen = floor(9000);

%% Code

addpath(genpath(s_path));
addpath(genpath(Toolbox_path));

filesM1L = dir([s_path,'/M1L']);
filesS1L =  dir([s_path,'/S1L']);
filesS1R =  dir([s_path,'/S1R']);

filesM1L(1:2) = [];
filesS1L(1:2) = []; 
filesS1R(1:3) = [];

prompt = {'Choose graph'};
dlg_title = 'Information required';
num_lines = 1;
def = {'1'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
s_which = str2num(answer{1});


figure(2)
total = [];

if s_which == 1 || s_which == 5

    for index = 1:length(filesM1L) 
        info = openSpike2(filesM1L(index).name);
        [p,q] = rat(info(whereisM1L).header.sampleinterval/s_SI,0.001);
        info(whereisM1L).data = resample(double(info(whereisM1L).data),p,q);
        total = [total; info(whereisM1L).data(1:s_minLen)'];
        time = 0:0.0001*s_SI:0.0001*s_SI*(length(info(whereisM1L).data)-1);
        plot(time,info(whereisM1L).data,'r');
        xlabel('Temps (ms)');
        ylabel('Potentiel (uV)');
        hold on;
        
    end
end

if s_which == 2 || s_which == 4 || s_which == 5

    
    for index = 1:length(filesS1L) 
        info = openSpike2(filesS1L(index).name);
        [p,q] = rat(info(whereisS1L).header.sampleinterval/s_SI,0.001);
        info(whereisS1L).data = resample(double(info(whereisS1L).data),p,q);
        total = [total; info(whereisS1L).data(1:s_minLen)'];
        time = 0:0.0001*s_SI:0.0001*s_SI*(length(info(whereisS1L).data)-1);
        plot(time,info(whereisS1L).data,'r');
        xlabel('Temps (ms)');
        ylabel('Potentiel (uV)');
        hold on; 
    end
end
    
if s_which == 3 || s_which == 4 || s_which == 5

    
    for index = 1:length(filesS1R) 
        info = openSpike2(filesS1R(index).name);
        [p,q] = rat(info(whereisS1R).header.sampleinterval/s_SI,0.001);
        info(whereisS1R).data = resample(double(info(whereisS1R).data),p,q);
        total = [total; info(whereisS1R).data(1:s_minLen)'];
        time = 0:0.0001*s_SI:0.0001*s_SI*(length(info(whereisS1R).data)-1);
        plot(time,info(whereisS1R).data,'r');
        xlabel('Temps (ms)');
        ylabel('Potentiel (uV)');
        hold on;
        
    end
end

s_mean = mean(total,1);
time = 0:0.0001*s_SI:0.0001*s_SI*(length(s_mean)-1);
figure(2)
plot(time,s_mean,'k','LineWidth',2)
switch s_which
    case 1 
        title('M1L')
    case 2 
        title('S1L')
    case 3
        title ('S1R')
    case 4
        title ('S1L & S1R')
    case 5
        title ('S1L & S1R & M1R')
end