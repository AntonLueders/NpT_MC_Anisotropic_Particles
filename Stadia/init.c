#include "init.h"
#include "distance.h"
#include "verletlist.h"

// ----------------------------------------------------------------------------------------

// Initilizes the simulation
void InitRandom(Particle *P, double L) {    

    long counter_overlap = 0;
    
    for (int i = 0; i < N; i++) {   
		
		if(counter_overlap > (long)N * (long)N) {
            printf("ERROR: System too dense.\n");
            exit(1);
        }
		
		for(int d = 0; d < dim; d++) {
            P[i].r_old[d] = 0.;
            P[i].e_old[d] = 0.;
        }
		P[i].q = q;
		P[i].sigma = sigma;
		
        for (int d = 0; d < dim; d++) {
                
            P[i].r[d] = L * gsl_rng_uniform(generator);
		}
		
		double rand_angle = 2. * M_PI * gsl_rng_uniform(generator);
		P[i].e[0] = cos(rand_angle);
		P[i].e[1] = sin(rand_angle);

		for (int j = 0; j < i; j++) {
            
			double disij = DistanceOverlap(&(P[i]), &(P[j]), L);     // See distance.c
			double disji = L;
			if (disij <= 3. / 2. * sigma) {
				disji = DistanceOverlap(&(P[j]), &(P[i]), L);	// See distance.c
			}
            if(disij < sigma || disji < sigma) {
                i--;
                counter_overlap++;
                break;
            }
        }
	}
	
	for (int i = 0; i < N; i++) {
		P[i].verList = malloc(sizeof(Node));
        P[i].verList->value = i;
        P[i].verList->next = NULL;	
	}
	
	printf("Init. successful.\n");
}

// ----------------------------------------------------------------------------------------
