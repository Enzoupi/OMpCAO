clear all
close all

!lecture du fichier
M=dlmread('Nbthreads.dat')

!Plot
figure(1)
hold on
plot(M(:,1),M(:,2),'-+')
plot(M(:,1),M(:,3),'-+')
plot(M(:,1),M(:,4),'-+')
plot(M(:,1),M(:,5),'-+')
xlabel('Nombre d éléments dans le vecteur à trier')
ylabel('Temps de calcul en secondes')
legend('1 Thread','2 Threads','3 Threads','4 Threads')
title('Temps de calcul en secondes en fonction du nombre de threads et de la taille du vecteur à trier')