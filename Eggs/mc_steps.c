#include "mc_steps.h"
#include "distance.h"

// ----------------------------------------------------------------------------------------

// Translational random step: It creats a test particles with the properties of
// the chosen particle. Then it displace this test particle randomly. 
// If the test particle does not overlap with another one, the step is accepted.
int MetropolisTranslation(int index, Particle *P, double L) {
	
	Particle test_particle = P[index];
    	
	for (int d = 0; d < dim; d++) {
		test_particle.r[d] = PeriodicPos(test_particle.r[d] + delta * (gsl_rng_uniform(generator) - 0.5), L);
	}
	
	Node *current = test_particle.verList->next;
    while (current != NULL) {
		
		int i = current->value;
        
		if (i != index) {
			double disij = DistanceOverlap(&test_particle, &(P[i]), L);     // See distance.c
			double disji = 0.;
			if (disij < b) {
				disji = Distance(&(P[i]), &test_particle, L, 0);     // See distance.c
			}
			if(disij < 0. || disji < 0.) {
								
				return 0;
			}
		}
		current = current->next;
    } 
	
	for (int d = 0; d < dim; d++) {
		P[index].r[d] = test_particle.r[d];
	}	
	return 1;
}

// ----------------------------------------------------------------------------------------

// Rotational random step: It creats a test particles with the properties of
// the chosen particle. Then it rotates this test particle randomly. 
// If the test particle does not overlap with another one, the step is accepted.
int MetropolisRotation(int index, Particle *P, double L) {
	
	double rand_unity_vector[dim];

        double rand_angle = 2. * M_PI * gsl_rng_uniform(generator);
        rand_unity_vector[0] = cos(rand_angle);
        rand_unity_vector[1] = sin(rand_angle);

        Particle test_particle = P[index];

        double esq = 0.;
        for (int d = 0; d < dim; d++) {
                test_particle.e[d] = test_particle.e[d] + delta_r * rand_unity_vector[d];
                esq += test_particle.e[d] * test_particle.e[d];
        }

        double esq_norm = sqrt(esq);

        for (int d = 0; d < dim; d++) {
                test_particle.e[d] /= esq_norm;
        }

	Node *current = test_particle.verList->next;
    while (current != NULL) {
		
		int i = current->value;
            
		if (i != index) {
			double disij = DistanceOverlap(&test_particle, &(P[i]), L);     // See distance.c
			double disji = 0.;
			if (disij < b) {
				disji = Distance(&(P[i]), &test_particle, L, 0); 
			}
			if(disij < 0. || disji < 0.) {
				
				return 0;
			}
		}
		current = current->next;
    }
	
	for (int d = 0; d < dim; d++) {
		P[index].e[d] = test_particle.e[d];
	}
	
	return 1;
}

// ----------------------------------------------------------------------------------------

// Random volume step for the implementation of the NpT ensemble
double VolumeMCStep(double L, Particle *P) {
	
	double V = pow(L, (double)dim);
	double V_new = V + (gsl_rng_uniform(generator) - 0.5) * delta_V;
	double L_temp = pow(V_new, 1./ (double)dim);
	
	Particle test_particles[N];
	
	for (int i = 0; i < N; i++) {
		
		test_particles[i] = P[i];
		for (int d = 0; d < dim; d++) {
			
			test_particles[i].r[d] *= L_temp / L;
		}
		
		Node *current = test_particles[i].verList->next;
		while (current != NULL) {
				
			int j = current->value;
				
			if (j < i) {
				
				double disij = DistanceOverlap(&(test_particles[i]), &(test_particles[j]), L_temp);     // See distance.c
				double disji = 0.;
				if (disij < b) {				
					disji = Distance(&(test_particles[j]), &(test_particles[i]), L_temp, 0);     // See distance.c
				}
				if(disij < 0. || disji < 0.) {
					return L;
				}
			}
			current = current->next;
		}  
	}
	
	double p_acc = exp(- p / T * (V_new - V) + (double)N * log(V_new / V));
	if (1. < p_acc) {
		p_acc = 1.;
	}
	
	if (gsl_rng_uniform(generator) < p_acc) {
		
		for (int i = 0; i < N; i++) {
			for (int d = 0; d < dim; d++) {
				P[i].r[d] = test_particles[i].r[d];
			}			
		}
	
		return L_temp;
	} else {
		
		return L; 
	}
}

// ----------------------------------------------------------------------------------------
