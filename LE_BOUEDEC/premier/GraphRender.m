clear all
close all

M = dlmread('/homedir/UR1/13007651/Documents/M2/ProgPar/TP1_Rendu/premier/Nbthreads.dat')
figure(1)
hold on
plot(M(:,1),M(:,2),'-+')
plot(M(:,1),M(:,3),'-+')
title('Temps de calcul en fonction du nombre de threads et du scheduling')
legend('Static','Dynamic')
ylabel('Temps de calcul en secondes')
xlabel('Nombre de threads')
hold off


M2 = dlmread('/homedir/UR1/13007651/Documents/M2/ProgPar/TP1_Rendu/premier/Scheduling.dat')
figure(2)
hold on
plot(M2(:,1),M2(:,3),'-+')
plot(M2(:,1),M2(:,2),'-+')
title('Temps de calcul en fonction du la taille des chunks')
legend('Schedule = Static','Schedule=Dynamic')
ylabel('Temps de calcul en secondes')
xlabel('Taille des chunks')

