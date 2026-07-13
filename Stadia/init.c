#include "init.h"
#include "distance.h"
#include "verletlist.h"

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
		
		double esq = 0.;
		
        for (int d = 0; d < dim; d++) {
                
            P[i].r[d] = L * gsl_rng_uniform(generator);
			P[i].e[d] = 2. * (gsl_rng_uniform(generator) - 0.5);
			esq += P[i].e[d] * P[i].e[d];
		}
		
		for (int d = 0; d < dim; d++) {
			P[i].e[d] /= sqrt(esq);
		}

		for (int j = 0; j < i; j++) {
            
			double disij = DistanceOverlap(&(P[i]), &(P[j]), L);     // See distance.c
			double disji = L;
			if (disij <= 3. / 2. * sigma) {
				disji = DistanceOverlap(&(P[j]), &(P[i]), L);
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
