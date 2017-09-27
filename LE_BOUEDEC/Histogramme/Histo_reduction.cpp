#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>
#ifdef _OPENMP
#include "omp.h"
#endif

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

void FillTable(double *table, int N)
{
	//#pragma omp parallel for schedule(dynamic)
	for (int i=0 ; i<N ; i++)
	{
		table[i] = fRand(0.0,10.0);
		//std::cout << "table[" << i << "] = " << table[i] << std::endl;
	}	
}

void MakeHisto(int *histo, double *table, int N)
{
	#pragma omp parallel
	{
	int histo_private[10]={0};
	#pragma omp for schedule(dynamic)
	for (int i=0; i<N ; i++)
	{
		bool test = false;
		int j=0;
		while (test == false)
		{
			if (table[i]-j<1.0)
			{
				test = true;
				histo_private[j]++;
			}
			j++;
		}

	}
	for (int k=0; k<10;k++)
		{
			#pragma omp atomic
			histo[k]=histo[k] + histo_private[k];
		}
	}
	
}

void CheckHisto(int *histo)
{
	int sum = 0;
	for (int i=0; i<10 ; i++)
	{
		sum += histo[i];
		printf("histo[%d]=%d\n",i,histo[i]);
	}
	printf("somme elements dans histo : %d\n",sum);
}

int main ()
{
	// Random number seed initialisation
	srand(time(NULL));   // should only be called once
	
	std::cout << "Enter the size of table :" << std::endl;
	int N;
	std::cin >> N;
	double *table;
	table = new double [N];
	int histo[10]={0};
	FillTable(table, N);
	// get max number of threads
	int maxthreads = omp_get_max_threads();
	for (int nbcpu=1 ; nbcpu <= maxthreads ; nbcpu++)
	{
		omp_set_num_threads(nbcpu);
		double t0 = omp_get_wtime();
		MakeHisto(histo,table,N);
		double t1 = omp_get_wtime();
		CheckHisto(histo);
		std::cout << "nb cores = " << nbcpu << std::endl;
		for (int i=0; i<10 ; i++)
		{
			//std::cout << i << " apparait " << histo[i] << std::endl;
			histo[i]=0;
		}
		std::cout << "Time spent = " << t1-t0 << std::endl;
	}
	
	
}
