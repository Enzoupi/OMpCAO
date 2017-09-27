#include <unistd.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#ifdef _OPENMP
#include "omp.h"
#endif

int cmbpremier_s(int n)
{
    int nbpremier = 2;
#pragma omp parallel for schedule(static) reduction(+ : nbpremier)
    for (int i = 4; i < n; i++)
    {
        int j = 2;
        bool test = false;
        while (j < i and test == false)
        {
            if (i % j == 0)
            {
                test = true;
            }
            j += 1;
        }

        if (j == (i))
        {
            nbpremier += 1;
        }
    }

    return nbpremier;
}

int cmbpremier_s(int n, int chunksize)
{
    int nbpremier = 2;
#pragma omp parallel for schedule(static, chunksize) reduction(+ : nbpremier)
    for (int i = 4; i < n; i++)
    {
        int j = 2;
        bool test = false;
        while (j < i and test == false)
        {
            if (i % j == 0)
            {
                test = true;
            }
            j += 1;
        }

        if (j == (i))
        {
            nbpremier += 1;
        }
    }

    return nbpremier;
}

int cmbpremier_d(int n)
{
    int nbpremier = 2;
#pragma omp parallel for schedule(dynamic) reduction(+ : nbpremier)
    for (int i = 4; i < n; i++)
    {
        int j = 2;
        bool test = false;
        while (j < i and test == false)
        {
            if (i % j == 0)
            {
                test = true;
            }
            j += 1;
        }

        if (j == (i))
        {
            nbpremier += 1;
        }
    }

    return nbpremier;
}

int cmbpremier_d(int n, int chunksize)
{
    int nbpremier = 2;
#pragma omp parallel for schedule(dynamic, chunksize) reduction(+ : nbpremier)
    for (int i = 4; i < n; i++)
    {
        int j = 2;
        bool test = false;
        while (j < i and test == false)
        {
            if (i % j == 0)
            {
                test = true;
            }
            j += 1;
        }

        if (j == (i))
        {
            nbpremier += 1;
        }
    }

    return nbpremier;
}

int main()
{
    // get max number of threads
    int maxthreads = omp_get_max_threads();

    //prepare outputs
    FILE *fo1;
    fo1 = fopen("Nbthreads.dat", "w");
    FILE *fo2;
    fo2 = fopen("Scheduling.dat", "w");

    // Which tests to perform
    bool TEST_NB_THREADS = true;
    bool TEST_SCHEDULE = true;

    // How many integers to check
    int parameter = 50000;

    if (TEST_NB_THREADS)
    {
        for (int nbcpu = 1; nbcpu <= maxthreads; nbcpu++)
        {
            omp_set_num_threads(nbcpu);
            std::cout << "number of threads : " << nbcpu << std::endl;
            //time static
            int res_static;
            double t0 = omp_get_wtime();
            res_static = cmbpremier_s(parameter);
            double t1 = omp_get_wtime();
            double time_static = t1 - t0;
            //time_dynamic
            int res_dynamic;
            double t2 = omp_get_wtime();
            res_dynamic = cmbpremier_d(parameter);
            double t3 = omp_get_wtime();
            double time_dynamic = t3 - t2;
            std::cout << "time spent static :" << time_static << std::endl;
            std::cout << "result schedule static :" << res_static << std::endl;
            std::cout << "time spent dynamic :" << time_dynamic << std::endl;
            std::cout << "result schedule dynamic :" << res_dynamic << std::endl;
            fprintf(fo1, "%d %f %f\n", nbcpu, time_static, time_dynamic);
        }
    }

    if (TEST_SCHEDULE)
    {
        omp_set_num_threads(maxthreads);
        for (int chunksize = 1; 50; chunksize = chunksize + 1)
        {
            //static
            double t0 = omp_get_wtime();
            cmbpremier_s(parameter, chunksize);
            double t1 = omp_get_wtime();
            double time_static = t1 - t0;
            //dynamic
            double t2 = omp_get_wtime();
            cmbpremier_d(parameter, chunksize);
            double t3 = omp_get_wtime();
            double time_dynamic = t3 - t2;
            fprintf(fo2, "%d %f %f\n", chunksize, time_static, time_dynamic);
        }
    }

    return 0;
}
