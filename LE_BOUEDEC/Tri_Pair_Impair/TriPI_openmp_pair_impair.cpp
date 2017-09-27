#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef _OPENMP
#include "omp.h"
#endif

void TriPI(float *table, int dimension)
/* Fonction qui a pour but de trier le vecteur "table" qui a "dimension"
 * composantes. 
 * <-- dimension : dimension du vecteur à trier
 * <-> table : vecteur original à trier qui sera modifié */
{
    if (dimension % 2 != 0)
    {
        float temp;
        for (int loopcounter = 0; loopcounter < dimension; loopcounter++)
        {
            if (loopcounter % 2 == 0) //phase paire donc indices impairs
            {
#pragma omp parallel for schedule(static)
                for (int i = 1; i <= dimension-2; i = i + 2)
                {
                    if ((*(table + i)) > *(table + i + 1))
                    {
                        temp = *(table + i);
                        *(table + i) = *(table + i + 1);
                        *(table + i + 1) = temp;
                    }
                }
            }
            else //phase impaire donc indices pairs ici
            {
#pragma omp parallel for schedule(static)
                for (int i = 0; i < dimension-2; i = i + 2)
                {
                    if (*(table + i) > *(table + i + 1))
                    {
                        temp = *(table + i);
                        *(table + i) = *(table + i + 1);
                        *(table + i + 1) = temp;
                    }
                }
            }
        }
    }
    else
    {
        int loopcounter = 0;
        float temp;
        for (loopcounter = 1; loopcounter <= dimension; loopcounter++)
        {
            if (loopcounter % 2 == 0) //phase paire donc indices impairs
            {
#pragma omp parallel for schedule(static)
                for (int i = 1; i < dimension - 2; i = i + 2)
                {
                    if ((*(table + i)) > *(table + i + 1))
                    {
                        temp = *(table + i);
                        *(table + i) = *(table + i + 1);
                        *(table + i + 1) = temp;
                    }
                }
            }
            else //phase impaire donc indices pairs ici
            {
#pragma omp parallel for schedule(static)
                for (int i = 0; i < dimension; i = i + 2)
                {
                    if (*(table + i) > *(table + i + 1))
                    {
                        temp = *(table + i);
                        *(table + i) = *(table + i + 1);
                        *(table + i + 1) = temp;
                    }
                }
            }
        }
    }
}

void FillRandomVector(float *vec, int dims, float max)
{
#pragma omp parallel for
    for (int i = 0; i < dims; i++)
    {
        *(vec + i) = rand() * max;
    }
}

int main(void)
{
    // Random number seed initialisation
    srand(time(NULL)); // should only be called once

    //Which actions to perform
    bool TEST = true;
    bool NB_THREADS = true;

    // get max number of threads
    int maxthreads = omp_get_max_threads();

    //prepare outputs
    FILE *foutput;
    foutput = fopen("Nbthreads.dat", "w");

    if (TEST)
    {
        float table_even[10] = {5.3, 8.6, 2.4, 4.5, 1.2, -8.2, 7.6, 2.3, 4.1, -10.0};
        std::cout << "Original Test Table Even:" << std::endl;
        for (int i = 0; i < 10; i++)
        {
            std::cout << *(table_even + i) << std::endl;
        }
        TriPI(table_even, 10);
        std::cout << "Ordered Test Table :" << std::endl;
        for (int i = 0; i < 10; i++)
        {
            std::cout << *(table_even + i) << std::endl;
        }

        float table_odd[9] = {5.3, 8.6, 2.4, 4.5, 1.2, -8.2, 7.6, 2.3, -10.0};
        std::cout << "Original Test Table Odd:" << std::endl;
        for (int i = 0; i < 9; i++)
        {
            std::cout << *(table_odd + i) << std::endl;
        }
        TriPI(table_odd, 9);
        std::cout << "Ordered Test Table Odd :" << std::endl;
        for (int i = 0; i < 9; i++)
        {
            std::cout << *(table_odd + i) << std::endl;
        }
    }

    if (NB_THREADS)
    {
        int dimension_max = 50000;
        for (int k = 1000; k < dimension_max; k = k + 5001)
        {
            float *vec = new float[k];
            float max = 50.0;
            //Fill table of size "dimension"
            FillRandomVector(vec, k, max);
            fprintf(foutput, "%d ", k);
            for (int i = 0; i < maxthreads; i++)
            {
                omp_set_num_threads(i + 1);
                double t0 = omp_get_wtime();
                TriPI(vec, k);
                double t1 = omp_get_wtime();
                double time_spent = t1 - t0;
                std::cout << "number of threads :" << i + 1 << " time elapsed :" << time_spent << std::endl;
                fprintf(foutput, "%f ", time_spent);
            }
            fprintf(foutput, "\n");
        }
    }

    return 0;
}
