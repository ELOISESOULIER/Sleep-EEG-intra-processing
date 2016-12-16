%st_info : a file opened with function openSpike2 with channel 1 containing
%intracellular signal and channel 2 EEG signal
%m_upstates : 2-columns cell containing all the different upstate parts of
%the intra in one column and the corresponding part of the EEG in the other
%v_number : a vector containing the begining index of each upstate
%st_newInfo = info with subsampled vm

function [m_upstates,v_number,st_newInfo] = upstate_analysis(st_info)

%% data and parameters %%
s_alpha=2;%modulates the position of the threshold between up and down
s_mindistupdown = 5; %minimum acceptable distance between up and down state
st_newInfo=st_info;
v_EEG = st_info(2).data;   
v_Vm_origin  = st_info(1).data;
s_bord=30;

%setting intra and EEG to the same sample rate
s_subsampRate = st_info(2).header.sampleinterval/st_info(1).header.sampleinterval;
v_index = 1:s_subsampRate:length(v_Vm_origin);
st_newInfo(1).header.Sampling = st_info(1).header.Sampling/s_subsampRate;
st_newInfo(1).header.sampleinterval = st_info(1).header.sampleinterval*s_subsampRate;
s_sampInt = st_newInfo(1).header.sampleinterval;
v_Vm = v_Vm_origin(v_index);
s_min = min(length(v_EEG),length(v_Vm));
v_Vm = v_Vm(1:s_min); v_EEG = v_EEG(1:s_min);
st_newInfo(1).data = v_Vm;
st_newInfo(2).data = v_EEG;
assert(length(v_Vm)==length(v_EEG));


%% separate upstate and downstate %%

%figure(3)
%histogram(v_Vm,'BinWidth',0.5);
%title('distribution des valeurs de lintra');

%find indices of the 2 maxima assimilated to the means of the gaussian
%distributions of up and down values
[v_N,v_edges]=histcounts(v_Vm);
[~, s_index1] = max(v_N);
s_firstMax = v_edges(s_index1);
v_N(s_index1)=0;
[~, s_index2] = max(v_N);
s_secondMax = v_edges(s_index2);
%check that second max is not on the same gaussian as first max
while abs(s_secondMax - s_firstMax) < s_mindistupdown
    v_N(s_index2)=0;
    [~, s_index2] = max(v_N);
    s_secondMax = v_edges(s_index2);
end
s_downstate = min(s_firstMax,s_secondMax);
s_upstate =  max(s_firstMax,s_secondMax);
s_threshold = (s_upstate - s_downstate)/s_alpha;
%label up as "1" and down as "0"
v_label = v_Vm > (s_upstate - s_threshold);

%% Remove "bad" upstates %%

%too short
v_variation =  v_label(2: end) - v_label(1: end-1);
v_switch = find(v_variation);
for i=2:length(v_switch)
    if (v_switch(i) - v_switch(i-1)) < ceil(50/s_sampInt*1000) 
        v_label(v_switch(i-1): v_switch(i)) =...
            v_label(v_switch(i-1))*ones(1,v_switch(i)-v_switch(i-1)+1);
    end
end

%downstate < 100ms means next upstate is actually the end of the one
%before, and it should not be take  into account
v_variation =  v_label(2:end) - v_label(1:end-1);
v_switch = find(v_variation);
for i=2:length(v_switch)-1
    if v_variation(v_switch(i))==-1
        continue
        else
        if v_switch(i) - v_switch(i-1) < ceil(100/s_sampInt*1000)
            v_label(v_switch(i): v_switch(i+1)) =...
                zeros(1,v_switch(i+1)-v_switch(i)+1);
        end
    end    
end


%% Put all the upstates in cell m_upstates %%

%first approximation of the number of upstates
s_nbUpstates = sum(v_variation < -0.5) + v_label(1) + 1;

v_number = [];
m_upstates = cell(s_nbUpstates,2);
s_count = 0;
for i = 1:length(v_label)
   if v_label(i) == 1
        if i == 1
            s_count = 1;
            v_number = [1];
        else
        s_count = s_count + (v_label(i-1) == 0) ;
            if (v_label(i-1) == 0)
                v_number = [v_number,i];
            end
        end
        m_upstates{s_count,1} = [m_upstates{s_count,1},v_Vm(i)];
        m_upstates{s_count,2} = [m_upstates{s_count,2},v_EEG(i)];
   end 
end
%restrict m_upstates to the actual number of upstates
m_upstates=m_upstates(1:numel(v_number),:);

%% precise re-cutting of each upstate%%

for s_which = 1:numel(v_number)
    v_origin = m_upstates{s_which,1};
    if length(v_origin) > 100 %will never be used otherwise so why bother
        %low-pass filter
        v_auxG = 2*v_origin(1) - v_origin(end:-1:1);
        v_auxR = 2*v_origin(end)-v_origin(end:-1:1);
        m_up=[v_auxG,v_origin,v_auxR];
        v_passe_bas	= f_DesignIIRfilter(st_newInfo(1).header.Sampling,70,71);
        m_upbis = f_FilterIIR(m_up,v_passe_bas);
        m_upbis=m_upbis(floor(end/3)+1:floor(2*end/3));
        %cut at first and last local max
        v_diff= m_upbis(2:end) - m_upbis(1:end-1);
        A=find(v_diff<0);
        B=find(v_diff>0);
        %check that none of the maxima is noise
        inf=1;
        while m_upbis(A(inf))< s_upstate-(s_upstate - s_downstate)/5 && inf<length(A);
           inf=inf+1;
        end
        sup = length(B);
        while m_upbis(B(sup))< s_upstate- (s_upstate - s_downstate)/5 && sup>1;
           sup=sup-1;
        end
        m_upstates{s_which,1} = v_origin(A(inf):B(sup));
        m_upstates{s_which,2} = m_upstates{s_which,2}(A(inf):B(sup));
        v_number(s_which) = v_number(s_which)+A(inf);
        %fprintf('perte : %1f\n',(length(m_origin) -  length(m_origin(A(1):B(end))))/length(m_origin)*100);
        %affichage de la distribution avec d�coupage
        %figure(45)
        %[v_N,v_edges]=histcounts(v_Vm);
        %v_N = v_N/length(v_Vm);
        %bar(v_edges(1:end-1),v_N);
        %hold on;
        %[v_N,v_edges]=histcounts(m_upstates{s_which,1});
        %v_N = v_N/2/length(m_upstates{s_which,1});
        %bar(v_edges(1:end-1),v_N,'r');
        %hold off;
    else continue
    end
   
end

%% smoothen the edge-cutting %%
% for s_which=1:numel(v_number)
%     s_len = length(m_upstates{s_which,1});
%     v_leftvm = v_Vm(max(1,v_number(s_which)- s_bord):v_number(s_which));
%     v_rightvm = v_Vm(v_number(s_which)+s_len: min(end,v_number(s_which)+s_len+s_bord));
%     m_upstates{s_which,1} = [v_leftvm,m_upstates{s_which,1},v_rightvm];
%     v_leftEEG = v_EEG(max(1,v_number(s_which)- s_bord):v_number(s_which));
%     v_rightEEG = v_EEG(v_number(s_which)+s_len: min(end,v_number(s_which)+s_len+s_bord));
%     m_upstates{s_which,2} = [v_leftEEG,m_upstates{s_which,2},v_rightEEG];
%     v_number(s_which) = v_number(s_which) - s_bord;
% end

%% visualize %%
% figure(22)
% subplot(2,1,1)
% plot(v_Vm)
% hold on;
% for i = 1:length(v_number)
%     plot(v_number(i):v_number(i)+length(m_upstates{i,1})-1,m_upstates{i,1})
%     hold on
% end
% hold off;
% subplot(2,1,2)
% plot(v_EEG)
% 
% keyboard;
end










