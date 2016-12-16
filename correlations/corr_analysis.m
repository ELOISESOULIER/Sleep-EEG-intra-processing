figure(1)
subplot(2,2,1);
histogram(correlations.corr_sp, 50);
%hold on;
%histogram(correlationsenveloppes.corr_sp,20);
%hold off;
%legend('signal brut', 'enveloppes')
title('Correlations Spindle');
xlabel('Coefficient de correlation');

subplot(2,2,2);
histogram(correlations.corr_gamLo,50);
%hold on;
%histogram(correlationsenveloppes.corr_gamLo,20);
%hold off;
%legend('signal brut', 'enveloppes')
title('Correlations Gamma 30-60');
xlabel('Coefficient de correlation');

subplot(2,2,3);
histogram(correlations.corr_gamMid,50);
%hold on;
%histogram(correlationsenveloppes.corr_gamMid,20);
%hold off;
%legend('signal brut', 'enveloppes')
title('Correlations Gamma 60-80');
xlabel('Coefficient de correlation');

subplot(2,2,4);
histogram(correlations.corr_gamHi,50);
%hold on;
%histogram(correlationsenveloppes.corr_gamHi,20);
%hold off;
%legend('signal brut', 'enveloppes')
title('Correlations Gamma 80-100');
xlabel('Coefficient de correlation');


figure(2)
subplot(2,2,1);
histogram(correlations.timelag_sp,100);
title('Décalage Spindle');
subplot(2,2,2);
histogram(correlations.timelag_gamLo,100);
title('Décalage Gamma 30-60');
subplot(2,2,3);
histogram(correlations.timelag_gamMid,100);
title('Décalage Gamma 60-80');
subplot(2,2,4);
histogram(correlations.timelag_gamHi,100);
title('Décalage Gamma 80-100');

figure(3)
subplot(2,2,1);
plot(correlations.length,correlations.corr_sp,'*')
xlabel('longueur upstates, pts')
ylabel('correlation vm EEG')
title('Spindle');
subplot(2,2,2);
plot(correlations.length,correlations.corr_gamLo,'*')
xlabel('longueur upstates, pts')
ylabel('correlation vm EEG')
title('Gamma 30-60');
subplot(2,2,3);
plot(correlations.length,correlations.corr_gamMid,'*')
xlabel('longueur upstates, pts')
ylabel('correlation vm EEG')
title('Gamma 60-80');
subplot(2,2,4);
plot(correlations.length,correlations.corr_gamHi,'*')
xlabel('longueur upstates, pts')
ylabel('correlation vm EEG')
title('Gamma 80-100');

