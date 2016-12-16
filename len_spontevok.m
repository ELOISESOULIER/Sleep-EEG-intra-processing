%compute length of spontaneous and evoked upstates

v_filesEvok = [{'c976.4-fqcy.smr'};{'c1036.2_fqcy.smr'};{'c1046.5-fqcy.smr'};...
     {'c1054.2-fqcy.smr'};{'c1307-fqcy.smr'};{'c1520.9-fqcy.smr'};...
     {'c2051_fqcy.smr'};{'c1176.3-evokfqcy1.smr'};{'c1277.6-fqcy2Evok.smr'};...
     {'c1615.5-fqcyevok.smr'};{'c2098.4-fqcyEvok.smr'}];
 v_filesSpont = [{'c976.4-fqcy.smr'};{'c1036.2_fqcy.smr'};{'c1046.5-fqcy.smr'};...
     {'c1054.2-fqcy.smr'};{'c1307-fqcy.smr'};{'c1520.9-fqcy.smr'};...
     {'c2051_fqcy.smr'};{'c1176.3-spontfqcy.smr'};{'c1277.6-fqcySpont.smr'};...
     {'c1615.5-fqcyspont.smr'};{'c2098.4-fqcySpont.smr'}];
 
 assert(length(v_filesEvok) == length(v_filesSpont));
 s_nbFiles = length(v_filesEvok);
 
 v_lenSpont = [];
 v_lenEvok = [];
 
for file = 1:s_nbFiles
    
    info = openSpike2(v_filesEvok{file});
    [~,m_evok,~] = spont_evok_extract(info,false,true);
    info = openSpike2(v_filesSpont{file});
    [m_spont,~,info] = spont_evok_extract(info,true,false);
    
    for s_which = 1:length(m_spont)
        v_lenSpont = [v_lenSpont,length(m_spont{s_which,1})];
    end
    for s_which = 1:length(m_evok)
        v_lenEvok = [v_lenEvok,length(m_evok{s_which,1})];
    end
    
end

figure(10)
histogram(v_lenSpont,20);
hold on;
histogram(v_lenEvok,20);
hold off;
legend('Spontané','Évoqué')
