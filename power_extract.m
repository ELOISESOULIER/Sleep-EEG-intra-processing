% Extracts powere in predefinite frequency bands corresponding to two
% ranges of spindle and three of gamma frequency
%To extract any frequency band, use single_power_extract.m


function [v_spLoSpont,v_spLoEvok,v_spHiSpont,v_spHiEvok,v_gamLoSpont,v_gamLoEvok,...
    v_gamMidSpont,v_gamMidEvok,v_gamHiSpont,v_gamHiEvok] = ...
    power_extract(m_TFAbsVmSpont,m_TFAbsVmEvok,v_FreqAxis)

v_FreqAxis = v_FreqAxis(end:-1:1); 
m_TFAbsVmSpont = m_TFAbsVmSpont(end:-1:1,:);
m_TFAbsVmEvok = m_TFAbsVmEvok(end:-1:1,:);
index = 1;
v_spLoSpont = zeros(1,size(m_TFAbsVmSpont,2));
v_spHiSpont = zeros(1,size(m_TFAbsVmSpont,2));
v_gamLoSpont = zeros(1,size(m_TFAbsVmSpont,2));
v_gamMidSpont = zeros(1,size(m_TFAbsVmSpont,2));
v_gamHiSpont = zeros(1,size(m_TFAbsVmSpont,2));

v_spLoEvok = zeros(1,size(m_TFAbsVmEvok,2));
v_spHiEvok = zeros(1,size(m_TFAbsVmEvok,2));
v_gamLoEvok = zeros(1,size(m_TFAbsVmEvok,2));
v_gamMidEvok = zeros(1,size(m_TFAbsVmEvok,2));
v_gamHiEvok = zeros(1,size(m_TFAbsVmEvok,2));

while v_FreqAxis(index)<9
   index = index + 1;
end

start = index;
while v_FreqAxis(index)<12
    v_spLoSpont = v_spLoSpont + m_TFAbsVmSpont(index,:);
    v_spLoEvok = v_spLoEvok + m_TFAbsVmEvok(index,:);
    index = index + 1;
end
v_spLoSpont = v_spLoSpont/(index - start);
v_spLoEvok = v_spLoEvok/(index - start);

start = index;
while v_FreqAxis(index)<16
    v_spHiSpont = v_spHiSpont + m_TFAbsVmSpont(index,:);
    v_spHiEvok = v_spHiEvok + m_TFAbsVmEvok(index,:);
    index = index + 1;
end
   
v_spHiSpont = v_spHiSpont/(index - start);
v_spHiEvok = v_spHiEvok/(index - start);

while v_FreqAxis(index)<30
    index = index + 1;
end

start = index;
while v_FreqAxis(index)<60
    v_gamLoSpont = v_gamLoSpont + m_TFAbsVmSpont(index,:);
    v_gamLoEvok = v_gamLoEvok + m_TFAbsVmEvok(index,:);
    index = index + 1;
end
   
v_gamLoSpont = v_gamLoSpont/(index - start);
v_gamLoEvok = v_gamLoEvok/(index - start);


start = index;
while v_FreqAxis(index)<80
     v_gamMidSpont = v_gamMidSpont + m_TFAbsVmSpont(index,:);
    v_gamMidEvok = v_gamMidEvok + m_TFAbsVmEvok(index,:);
    index = index + 1;
end
   
v_gamMidSpont = v_gamMidSpont/(index - start);
v_gamMidEvok = v_gamMidEvok/(index - start);

start = index;
while v_FreqAxis(index)<100
     v_gamHiSpont = v_gamHiSpont + m_TFAbsVmSpont(index,:);
    v_gamHiEvok = v_gamHiEvok + m_TFAbsVmEvok(index,:);
    index = index + 1;
end
   
v_gamHiSpont = v_gamHiSpont/(index - start);
v_gamHiEvok = v_gamHiEvok/(index - start);
