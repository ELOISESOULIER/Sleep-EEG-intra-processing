% cette fonction permet d'�limier les potentiels d'action par contraction
% de la d�riv�e au niveau du potentiel d'action

%ne fonctionne que pour 1 potentiel, utiliser traitementpa_general pour
%plusieurs potentiels d'action


function v_signal_traite = traitementpa(v_signal,s_sampleinterval)

	s_delta = 6; %largeur de la deformation autour du pa		
	
	% Pour éviter les "fausses" détection de pa dues à des pentes raides en
	% dessous du seuil
	if (max(v_signal)<-30)
		fprintf('pas de potentiel daction');
		v_signal_traite = v_signal;
	else
        

	%détecter la zone à modifier
        v_diff = (v_signal(2:end)-v_signal(1:end-1))/s_sampleinterval;
        [m,i] = max(v_diff);
        s_coupeinf = i -s_delta;
        [n,j] = min(v_diff);
        s_coupesup = j+2*s_delta;

        %modifier la derivee : shrink normalisé
        v_diff_traitee = v_diff;
        s_amp_moy = (mean(v_diff(1:i))+mean(v_diff(j:end)));
        v_diff_traitee(s_coupeinf:s_coupesup) = s_amp_moy./(1+abs(v_diff(s_coupeinf:s_coupesup)));
        %v_diff_traitee(s_coupeinf:s_coupesup) = v_diff(s_coupeinf:s_coupesup).*((s_amp_moy./v_diff(s_coupeinf:s_coupesup)).^5);

        % 
        % figure(2)
        % Ax(1)=subplot(2,1,1);
        % plot(v_diff);
        % Ax(2)=subplot(2,1,2);
        % plot(v_diff_traitee);
        % linkaxes(Ax, 'x');
	% title('Dérivée traitée');


        %reconstruire le signal
        s_ldt = length(v_diff_traitee);
        v_signal_traite = [v_signal(1)];

        for i=1:s_ldt
            s_elt = v_signal_traite(i) + s_sampleinterval*v_diff_traitee(i);
            v_signal_traite = [v_signal_traite,s_elt];
        end
        for i = s_coupeinf:s_coupesup
            factor = (i - s_coupeinf)/(s_coupesup - s_coupeinf);
            v_signal_traite(i) = v_signal_traite(i) +...
                (v_signal(i+s_coupesup -s_coupeinf)- v_signal_traite(i+s_coupesup-s_coupeinf))*factor;%(s_ldt-i)*2/s_ldt;
        end
        for i = s_coupesup:s_ldt
            v_signal_traite(i+1) = v_signal_traite(i) + s_sampleinterval*v_diff_traitee(i);
        end

 
        %figure(3)
        %plot(v_signal);
        %hold on;
        %plot(v_signal_traite);
	%hold off;

        %figure(4)
        %subplot(1,2,1)
        % [m_pa, v_TimeAxis, v_FreqAxis] = f_GaborTransformWait(v_signal,1/0.36*1000,5,100);
        % v_Clim = [prctile(m_pa(:),5),prctile(m_pa(:),95)];
        %f_ImageArray(m_pa, v_TimeAxis, v_FreqAxis,'colormap','jet','limits',v_Clim);
        %subplot(1,2,2)
        % [m_npa, v_TimeAxis, v_FreqAxis] = f_GaborTransformWait(v_signal_traite,1/0.36*1000,5,100);
        %f_ImageArray(m_npa, v_TimeAxis, v_FreqAxis,'colormap','jet','limits',v_Clim);
	%title('Effet de la modification');
    end
end
