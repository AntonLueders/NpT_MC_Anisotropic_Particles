#include "distance.h"

// ----------------------------------------------------------------------------------------

// Periodic boundary conditions (for position)
inline double PeriodicPos(double dr, double L) {

    double dr_new = dr - L * floor(dr / L);

    return dr_new;
}

// ----------------------------------------------------------------------------------------

// Periodic boundary conditions (for distance)
inline double PeriodicDis(double dr, double L) {

    double dr_new = dr - L * nearbyint(dr / L);

    return dr_new;
}

// ----------------------------------------------------------------------------------------

// This is the first function which is used to check if the particles are even close enouth
// that checks for overlap are necessary.
double DistanceOverlap(Particle *p_i, Particle *p_j, double L) {
	
	double rijsq = 0.;
	for (int d = 0; d < dim; d++) {
		rijsq += PeriodicDis(p_j->r[d] - p_i->r[d], L) * PeriodicDis(p_j->r[d] - p_i->r[d], L);
	}
	
	double disij = 0.;
	if (sqrt(rijsq) < 1.05 * q * sigma) {
		disij = Distance(p_i, p_j, L);
	} else {
		disij = sqrt(rijsq);
	}
	
	return disij;
}

// ----------------------------------------------------------------------------------------

double Sign(double value, double sign) {
	
	if (fabs(sign) <= 1e-10) {
		return 0.;
	}
	return value * sign / fabs(sign);
}

// ----------------------------------------------------------------------------------------

// Distance between spherocylinders as implemented in the program written for
// SIEMS, Ullrich, 2017. Computersimulationen von Modell-Kolloiden in Kanälen [Dissertation]. Konstanz: Universität Konstanz
double Distance(Particle *p_i, Particle *p_j, double L) {
	
	double rijei = 0.;
	double rijej = 0.;
	double eiej = 0.;
	double l_half_i = (p_i->q - 1.) * p_i->sigma / 2.;
	double l_half_j = (p_j->q - 1.) * p_j->sigma / 2.;
	
	double rij[dim];
	
	for (int d = 0; d < dim; d++) {
		rij[d] = PeriodicDis(p_j->r[d] - p_i->r[d], L);
	}
	
	for (int d = 0; d < dim; d++) {
		
		rijei += rij[d] * p_i->e[d];
		rijej += rij[d] * p_j->e[d];
		eiej += p_i->e[d] * p_j->e[d];
	}
	
	double denom_temp = 1.0 - eiej * eiej;
	double lambda_prime;
	double mu_prime;
	if (denom_temp <= 1e-10) {
		lambda_prime = rijei / 2.;
		mu_prime = -rijej / 2.;
	} else {
		lambda_prime = (rijei -eiej * rijej) / denom_temp;
		mu_prime = (-rijej + eiej * rijei) / denom_temp;
	}
	if (fabs(lambda_prime) > l_half_i || fabs(mu_prime) > l_half_j) {
		double aux1 = fabs(lambda_prime) - l_half_i;
		double aux2 = fabs(mu_prime) - l_half_j;
		if (aux1 > aux2) {
			lambda_prime = Sign(l_half_i, lambda_prime);
			mu_prime = lambda_prime * eiej - rijej;
			if (fabs(mu_prime) > l_half_j) {
				mu_prime = Sign(l_half_j, mu_prime);
			}
		} else {
			mu_prime = Sign(l_half_j, mu_prime);
			lambda_prime = mu_prime * eiej + rijei;
			if(fabs(lambda_prime) > l_half_i) {
				lambda_prime = Sign(l_half_i, lambda_prime);
			}
		}
	}
	
	double dissq = 0.;
	double dij[dim];
	for (int d = 0; d < dim; d++) {
		dij[d] = rij[d] + mu_prime * p_j->e[d] - lambda_prime * p_i->e[d];
		dissq += dij[d] * dij[d];
	}
	
	return sqrt(dissq);
}

// ----------------------------------------------------------------------------------------
