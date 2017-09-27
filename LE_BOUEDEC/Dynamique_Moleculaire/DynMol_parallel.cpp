#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "omp.h"

//
// Computation parameters
//

const int numDimensions = 3;
const int numParticles = 2000;

const double t_max = 0.05;
const double dt = 0.0001;


// Interaction potential

void computePotential (const double x, double& v, double& dv) {
    double y = (x<M_PI_2) ? x : M_PI_2;
    v = sin(y)*sin(y);
    dv = sin(2.0*y);
}


//
// Particle class
//

class Particle {

public:

    Particle () {};
    ~Particle () {};

// Public methods

    void init (const double* p, const double* v, const double* a){
        for (int i=0; i<numDimensions; i++) {
            position[i] = p[i];
            velocity[i] = v[i];
            acceleration[i] = a[i];
        }
        kinetic_energy = 0.0;
        resetForce();
        mass = 1.0;
        };

    void update () {
        for (int i=0; i<numDimensions; i++) {
            position[i] += velocity[i]*dt + 0.5*acceleration[i]*dt*dt;
            velocity[i] += 0.5*dt*(force[i]/mass + acceleration[i]);
            acceleration[i] = force[i]/mass;
        }
    };

    void computeForce (const Particle& other) {
        double distance;
        double rij[numDimensions];
        computeSeparation(other,rij,distance);
        double v;
        double dv;
        computePotential(distance,v,dv);
        potential_energy += 0.5*v;
        for (int i=0; i<numDimensions; i++)
            force[i] -= rij[i]*dv/distance;
    };

    void resetForce () {
        potential_energy = 0.0;
        for (int i=0; i<numDimensions; i++)
            force[i] = 0.0;
    };

    double getKineticEnergy () {
        computeKineticEnergy();
        return kinetic_energy;
    };

    double getPotentialEnergy () {
        return potential_energy;
    };

    void print () const {
        printf("Position: [");
        for (int i=0; i<numDimensions; i++)
            printf("%8.5f",position[i]);
        printf("]    ");
        printf("Velocity: [");
        for (int i=0; i<numDimensions; i++)
            printf("%8.5f",velocity[i]);
        printf("]    ");
        printf("Acceleration: [");
        for (int i=0; i<numDimensions; i++)
            printf("%8.5f",acceleration[i]);
        printf("]\n");
    };

private:

// Private methods

    double computeKineticEnergy () {
        kinetic_energy = 0.0;
        for (int i=0; i<numDimensions; i++)
            kinetic_energy += velocity[i]*velocity[i];
        kinetic_energy *= 0.5*mass;
    };

    void computeSeparation (const Particle& other, double* rij, double& distance) {
        distance = 0.0;
        if (this == &other) {
            for (int i=0; i<numDimensions; i++)
                rij[i] = 0.0;
        } else {
            for (int i=0; i<numDimensions; i++) {
                double h = position[i] - other.position[i];
                rij[i] = h;
                distance += h*h;
            }
            distance = sqrt(distance);
        }
    };

// Private data

    double mass;
    double kinetic_energy;
    double potential_energy;
    double position[numDimensions];
    double velocity[numDimensions];
    double acceleration[numDimensions];
    double force[numDimensions];

};


//
// Main program
//

int main () {

//Prepare output file
FILE *foutput;
foutput = fopen("Nb_cores_VS_Wall_time.dat","w");

// Number of cores :
int NbMaxThreads = omp_get_max_threads();
for (int k=1 ; k<=NbMaxThreads ; k++)
{
omp_set_num_threads(k);

// Size of box containing particles

    double box[numDimensions];
    for (int i=0; i<numDimensions; i++)
        box[i] = 10.0;

// Initialize the positions of each particle to a random point
// Set initial velocities and accelerations to 0

    double zero[numDimensions];
    for (int i=0; i<numDimensions; i++)
        zero[i] = 0.0;

    Particle particle[numParticles];

    for (int p=0; p<numParticles; p++) {
        double pos[numDimensions];
        for (int i=0; i<numDimensions; i++)
            pos[i] = box[i] * (double)(p+1)/numParticles;
//            pos[i] = box[i] * (double)rand()/RAND_MAX;
        particle[p].init(pos,zero,zero);
    }

// Loop over time:

    double ke = 0.0;
    double pe = 0.0;

    double t = 0.0;

    int count = 0;

    double t0 = omp_get_wtime();
    while (t < t_max) {

// Compute forces acting on each particle and their energies

        ke = 0.0;
        pe = 0.0;

#pragma omp parallel for reduction(+:ke,pe) schedule(dynamic)
        for (int p=0; p<numParticles; p++) {
            particle[p].resetForce();
            for (int other=0; other<numParticles; other++) {
                if (other != p) {
                    particle[p].computeForce(particle[other]);
                }
            }
            ke += particle[p].getKineticEnergy();
            pe += particle[p].getPotentialEnergy();
        }

// Compute new positions and velocities of each particle
// from forces acting on them

        for (int p=0; p<numParticles; p++)
            particle[p].update();

// Print current energy values

        if (count%100 == 0)
            printf("Time t = %8.5f   pe = %14.5f   ke = %14.5f   total = %14.5f\n",
                    t,pe,ke,pe+ke);

// Update t

        t += dt;
        count ++;
    }
    double t1 = omp_get_wtime();

    printf("Time t = %8.5f   pe = %14.5f   ke = %14.5f   total = %14.5f\n",
            t,pe,ke,pe+ke);

    printf("Time taken was %8.5f seconds\n",t1-t0);
    
    fprintf(foutput,"%d %f\n",k,t1-t0);
    
}

    return 0;
}
