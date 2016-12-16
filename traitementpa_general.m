% cette fonction permet d'�limier les potentiels d'action par contraction
% de la d�riv�e au niveau du potentiel d'action, sur plusieurs potentiels d'action

% s_distmin = distance minimale sur le signal entre 2 pa

function v_signal_traite = traitementpa_general(v_signal,s_dist_min,s_sample_interval)

	s_delta = floor(0.8*s_dist_min/s_sample_interval); %largeur de la deformation autour du pa
	v_signal_traite = v_signal;

	v_diff = (v_signal(2:end) - v_signal(1:end-1))/s_sample_interval;
	v_Index = find(v_diff >10); %r�cup�rer les 
	%enlever les doublons (deux points pour un seul pa)
	v_aux = (v_Index(2:end) - v_Index(1:end-1))>2;
	v_Index = [v_Index(v_aux),v_Index(end)];

	for i=1:length(v_Index)
	    s_current = v_Index(i);
	    s_coupeinf = s_current -s_delta
	    s_coupesup = s_current+s_delta
	    
	    v_signal_traite(s_coupeinf:s_coupesup) = traitementpa(v_signal(s_coupeinf:s_coupesup),s_sample_interval*1000);
	    
    end

%% visualiser l'effet sur une carte temps fréquence
    
	%[m_GaborWTpa, v_TimeAxis, v_FreqAxis] =f_GaborTransformWait(v_signal,info(1).header.Sampling,5,120);
	%v_Clim = [prctile(m_GaborWTpa(:),5),prctile(m_GaborWTpa(:),95)];
	%figure(12)
	%f_ImageArray(m_GaborWTpa, v_TimeAxis, v_FreqAxis,'colormap','jet','limits',v_Clim)


	%[m_GaborWTpa2, v_TimeAxis, v_FreqAxis] =f_GaborTransformWait(info(1).data,info(1).header.Sampling,5,120);
	%figure(22)
	%f_ImageArray(m_GaborWTpa, v_TimeAxis, v_FreqAxis,'colormap','jet','limits',v_Clim)
end
