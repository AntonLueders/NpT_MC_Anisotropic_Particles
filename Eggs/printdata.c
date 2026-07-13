#include "printdata.h"

// ----------------------------------------------------------------------------------------

void PrintOvito(FILE *f, Particle *P, int step, double L) {
    
    fprintf(f, "%d\n", N);
    fprintf(f, "step = %d\n", step + 1);

	for (int i = 0; i < N; i++) {	
	
		double phi = atan2(P[i].e[1], P[i].e[0]);
	
		fprintf(f, "s1\t%f\t%f", P[i].r[0], P[i].r[1]);
		fprintf(f, "\t%f\t%f\t%f", phi / M_PI * 180., cos(phi / 2.), sin(phi / 2.));
		fprintf(f, "\t%f\t%f\t%f\t", P[i].a, P[i].b, P[i].b);
		fprintf(f, "%f\n", L);
	}
}

// ----------------------------------------------------------------------------------------