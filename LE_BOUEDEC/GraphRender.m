M = dlmread('Nbthreads.dat')

figure(1)
hold on
plot(M(:,1),M(:,2),'-+')
plot(M(:,1),M(:,3),'-+')
plot(M(:,1),M(:,4),'-+')
plot(M(:,1),M(:,5),'-+')
title('Temps de calcul en fonction du nombre de threads et de la taille la matrice')
legend('1 Thread','2 Threads','3 Threads','4 Threads')
ylabel('Temps de calcul en secondes')
xlabel('Taille des matrices carrées à multiplier')
hold off

figure(2)
hold on
M2 = dlmread('Schedule.dat')
plot(M2(:,1),M2(:,2),'-+')
plot(M2(:,1),M2(:,3),'-+')
title('Temps de calcul en fonction du la taille des chunks')
legend('Schedule = Static','Schedule = Dynamic')
ylabel('Temps de calcul en secondes')
xlabel('Taille des chunks')