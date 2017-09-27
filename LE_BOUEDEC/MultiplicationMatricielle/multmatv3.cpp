#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <stdlib.h>
#ifdef _OPENMP
#include "omp.h"
#endif

double RandMatMult(int n, int m, int p)
{
	clock_t begin = clock();
	double *aspace = new double[m * n];
	double **M1 = new double *[m];
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < m; i++)
	{
		M1[i] = &aspace[i * n];
	}

	double *bspace = new double[m * p];
	double **M2 = new double *[p];
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < p; i++)
	{
		M2[i] = &bspace[i * m];
	}

	double *cspace = new double[n * p];
	double **MR = new double *[p];
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < p; i++)
	{
		MR[i] = &cspace[i * n];
	}

#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			M1[i][j] = rand() % 15;
		}
	}

#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < p; j++)
		{
			M2[i][j] = rand() % 5;
		}
	}

	//Multiplication de matrices
	double t0 = omp_get_wtime();
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < p; j++)
		{
			for (int k = 0; k < m; k++)
			{
				MR[i][j] = MR[i][j] + M1[i][k] * M2[k][j];
			}
		}
	}
	double t1 = omp_get_wtime();
	double time_spent = t1 - t0;

	/*
	//Output pour verification
	for (int i=0 ; i<n ; i++)
	{
		for (int j=0 ; j<m ; j++)
		{
			printf(" %f ",M1[i][j]);
		}
	printf("\n");
	}

	for (int i=0 ; i<m ; i++)
	{
		for (int j=0 ; j<p ; j++)
		{
			printf(" %f ",M2[i][j]);
		}
	printf("\n");
	}

	for (int i=0 ; i<n ; i++)
	{
		for (int j=0 ; j<p ; j++)
		{
			printf(" %f ",MR[i][j]);
		}
	printf("\n");
	}
	*/

	return time_spent;
}

double RandMatMultDynamic(int n, int m, int p, int chunk)
{
	clock_t begin = clock();
	double *aspace = new double[m * n];
	double **M1 = new double *[m];
#pragma omp parallel for schedule(dynamic, chunk)
	for (int i = 0; i < m; i++)
	{
		M1[i] = &aspace[i * n];
	}

	double *bspace = new double[m * p];
	double **M2 = new double *[p];
#pragma omp parallel for schedule(dynamic, chunk)
	for (int i = 0; i < p; i++)
	{
		M2[i] = &bspace[i * m];
	}

	double *cspace = new double[n * p];
	double **MR = new double *[p];
#pragma omp parallel for schedule(dynamic, chunk)
	for (int i = 0; i < p; i++)
	{
		MR[i] = &cspace[i * n];
	}

#pragma omp parallel for schedule(dynamic, chunk)
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			M1[i][j] = rand() % 15;
		}
	}

#pragma omp parallel for schedule(dynamic, chunk)
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < p; j++)
		{
			M2[i][j] = rand() % 5;
		}
	}

	//Multiplication de matrices
	double t0 = omp_get_wtime();
#pragma omp parallel for schedule(dynamic, chunk)
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < p; j++)
		{
			for (int k = 0; k < m; k++)
			{
				MR[i][j] = MR[i][j] + M1[i][k] * M2[k][j];
			}
		}
	}
	double t1 = omp_get_wtime();
	double time_spent = t1 - t0;

	/*
	//Output pour verification
	for (int i=0 ; i<n ; i++)
	{
		for (int j=0 ; j<m ; j++)
		{
			printf(" %f ",M1[i][j]);
		}
	printf("\n");
	}

	for (int i=0 ; i<m ; i++)
	{
		for (int j=0 ; j<p ; j++)
		{
			printf(" %f ",M2[i][j]);
		}
	printf("\n");
	}

	for (int i=0 ; i<n ; i++)
	{
		for (int j=0 ; j<p ; j++)
		{
			printf(" %f ",MR[i][j]);
		}
	printf("\n");
	}
	*/

	return time_spent;
}

double RandMatMultStatic(int n, int m, int p, int chunk)
{
	clock_t begin = clock();
	double *aspace = new double[m * n];
	double **M1 = new double *[m];
#pragma omp parallel for schedule(static, chunk)
	for (int i = 0; i < m; i++)
	{
		M1[i] = &aspace[i * n];
	}

	double *bspace = new double[m * p];
	double **M2 = new double *[p];
#pragma omp parallel for schedule(static, chunk)
	for (int i = 0; i < p; i++)
	{
		M2[i] = &bspace[i * m];
	}

	double *cspace = new double[n * p];
	double **MR = new double *[p];
#pragma omp parallel for schedule(static, chunk)
	for (int i = 0; i < p; i++)
	{
		MR[i] = &cspace[i * n];
	}

#pragma omp parallel for schedule(static, chunk)
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			M1[i][j] = rand() % 15;
		}
	}

#pragma omp parallel for schedule(static, chunk)
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < p; j++)
		{
			M2[i][j] = rand() % 5;
		}
	}

	//Multiplication de matrices
	double t0 = omp_get_wtime();
#pragma omp parallel for schedule(static, chunk)
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < p; j++)
		{
			for (int k = 0; k < m; k++)
			{
				MR[i][j] = MR[i][j] + M1[i][k] * M2[k][j];
			}
		}
	}
	double t1 = omp_get_wtime();
	double time_spent = t1 - t0;

	/*
	//Output pour verification
	for (int i=0 ; i<n ; i++)
	{
		for (int j=0 ; j<m ; j++)
		{
			printf(" %f ",M1[i][j]);
		}
	printf("\n");
	}

	for (int i=0 ; i<m ; i++)
	{
		for (int j=0 ; j<p ; j++)
		{
			printf(" %f ",M2[i][j]);
		}
	printf("\n");
	}

	for (int i=0 ; i<n ; i++)
	{
		for (int j=0 ; j<p ; j++)
		{
			printf(" %f ",MR[i][j]);
		}
	printf("\n");
	}
	*/

	return time_spent;
}
int main()
{
	// Random number seed initialisation
	srand(time(NULL)); // should only be called once

	// get max number of threads
	int maxthreads = omp_get_max_threads();

	//prepare outputs
	FILE *foutput;
	foutput = fopen("Nbthreads.dat", "w");
	FILE *fo2;
	fo2 = fopen("Schedule.dat", "w");

	bool TEST_NB_THREADS = true;
	bool TEST_SCHEDULE = true;

	if (TEST_NB_THREADS)
	{
		printf("Test number of threads\n");
		//Loop on matrix size
		for (int param = 100; param < 905; param += 100)
		{
			printf("Matrix size : %d\n",param);
			fprintf(foutput, "%d ", param);
			//Loop on thread(s) number
			for (int nbcpu = 1; nbcpu <= maxthreads; nbcpu++)
			{
				omp_set_num_threads(nbcpu);
				int n = param;
				int m = param;
				int p = param;
				double time_spent = RandMatMult(n, m, p);
				fprintf(foutput, "%f ", time_spent);
			}
			fprintf(foutput, "\n");
		}
	}

	//Test on chunksize and schedule
	if (TEST_SCHEDULE)
	{
		printf("Test on scheduling\n");
		int sizemat = 500;
		omp_set_num_threads(2);
		for (int i = 1; i <= sizemat / 2; i = i + 10)
		{
			printf("chunk size = %d",i);
			fprintf(fo2, "%d ", i);
			double time_static = RandMatMultStatic(sizemat, sizemat, sizemat, i);
			double time_dynamic = RandMatMultDynamic(sizemat, sizemat, sizemat, i);
			fprintf(fo2, "%f %f\n", time_static, time_dynamic);
		}
	}

	return 0;
}
