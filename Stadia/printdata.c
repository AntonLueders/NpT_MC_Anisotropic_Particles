#include "printdata.h"

// ----------------------------------------------------------------------------------------

void PrintOvito(FILE *f, Particle *P, int step, double L) {
    
    fprintf(f, "%d\n", N);
    fprintf(f, "step = %d\n", step + 1);

	for (int i = 0; i < N; i++) {	
	
		double phi = atan2(P[i].e[1], P[i].e[0]);
		// Kombinierter Quaternion nach Rotation um 90 Grad um die Y-Achse und dann um 90 Grad um die Z-Achse
		double w =  0.5 * sqrt(1. + cos(phi));
		double x = sin(phi) / (4.*w);
		double y = -(cos(phi)+1.)/(4.*w);
		double z =  sin(phi)/(4.*w);

		fprintf(f, "s1\t%f\t%f", P[i].r[0], P[i].r[1]);
		fprintf(f, "\t%f\t%f\t%f\t%f", w,x,y,z);
		fprintf(f, "\t%f\t%f\t", P[i].sigma / 2., (P[i].q - 1.) * P[i].sigma);
		fprintf(f, "%f\n", L);
	}
}

// ----------------------------------------------------------------------------------------
