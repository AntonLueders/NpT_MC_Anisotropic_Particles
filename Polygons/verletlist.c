#include "verletlist.h"
#include "distance.h"

void GenerateVerletLists(Particle *P, double L, int step) {
	
	for (int i = 0; i < N; i++) {
		
		clearList(P[i].verList);
		
		for (int j = 0; j < i; j++) {
			
			double disij = DistanceOverlap(&(P[i]), &(P[j]), L);     // See distance.c
			
			double limit = (r_v + 0.5 / sin(M_PI / (double)m) * a);
			if (lambda > 1.) {
				limit *= lambda;
			}
			
			if (disij < 0.) {
				printf("ERROR: Unexpected overlap after step: %d.\n", step);
				exit(1);
			} else if(disij < limit) {
				push(P[i].verList, P[j].verList->value);
				push(P[j].verList, P[i].verList->value);					
			}
		}
		
		for (int d = 0; d < dim; d++) {
			
			P[i].r_old[d] = P[i].r[d];
			P[i].e_old[d] = P[i].e[d];
		}
	}
}

bool CheckVerlet(Particle *P, double L) {
	
	double dissq_max[2] = {0., 0.};
	double l_half_i = 0.5 / sin(M_PI / (double)m) * a;
	if (lambda > 1.) {
		l_half_i *= lambda;
	}
	
	for (int i = 0; i < N; i++) {
	
		double dissq = 0.;
		double dcross = 0.;
		double de;
		double dr;
		
		for (int d = 0; d < dim; d++) {
			de = l_half_i * (P[i].e[d] - P[i].e_old[d]);
			dr = PeriodicDis(P[i].r[d] - P[i].r_old[d], L);
			
			dissq += de * de + dr * dr;
			dcross += de * dr;
		}
		dissq += 2. * fabs(dcross);
		
		if (dissq > dissq_max[1]) {
			dissq_max[0] = dissq_max[1];
			dissq_max[1] = dissq;
		} else if (dissq > dissq_max[0]) {
			dissq_max[0] = dissq;
		}
	}
	
	return (sqrt(dissq_max[0]) + sqrt(dissq_max[1])) >= (r_v - delta);
}

void FreeVerLists(Particle *P) {

    for(int i = 0; i < N; i++) {
        clearList(P[i].verList);     // See lists.c
        free(P[i].verList);
    }
}
