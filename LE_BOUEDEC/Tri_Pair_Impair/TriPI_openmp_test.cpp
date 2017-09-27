#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef _OPENMP
#include "omp.h"
#endif

void TriPI(float *table, int dimension)
{
    if (dimension % 2 != 0)
    {
        "Le nombre d'éléments n'est pas pair !";
        exit(EXIT_FAILURE);
    }
    int loopcounter = 0;
    float temp;
    for (loopcounter = 1; loopcounter <= dimension; loopcounter++)
    {
#pragma omp parallel
{
        if (loopcounter % 2 == 0) //phase paire donc indices impairs
        {
#pragma omp for schedule(static) private(temp)
            for (int i = 1; i < dimension - 1; i = i + 2)
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
#pragma omp for schedule(static) private(temp)
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

int CheckOrderedTable(float *vec, int dims)
{
    for (int i=0;i<dims-1;i++)
    {
        if(*(vec+i) > *(vec +i +1))
        {
            printf("Table is not ordered at index %d\n",i);
            return 1;
        }
    }
    return 0;
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
        float tabletest[10] = {5.3, 8.6, 2.4, 4.5, 1.2, -8.2, 7.6, 2.3, 4.1, 8.2};
        std::cout << "Original Test Table :" << std::endl;
        for (int i = 0; i < 10; i++)
        {
            std::cout << *(tabletest + i) << std::endl;
        }
        TriPI(tabletest, 10);
        std::cout << "Ordered Test Table :" << std::endl;
        for (int i = 0; i < 10; i++)
        {
            std::cout << *(tabletest + i) << std::endl;
        }
    }

    if (NB_THREADS)
    {
        int dimension_max = 15000;
        for (int k = 10; k < dimension_max; k = k + 400)
        {
            float *vec = new float[k];
            float max = 50.0;
            //Fill table of size "dimension"
            FillRandomVector(vec, k, max);
            printf("Table is of size %d\n",k);
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
                CheckOrderedTable(vec,k);
            }
            fprintf(foutput, "\n");
            delete  [] vec;
            
        }
    }

    return 0;
}
