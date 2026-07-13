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

double DistanceOverlap(Particle *p_i, Particle *p_j, double L) {
	
	double rijsq = 0.;
	for (int d = 0; d < dim; d++) {
		rijsq += PeriodicDis(p_j->r[d] - p_i->r[d], L) * PeriodicDis(p_j->r[d] - p_i->r[d], L);
	}
	
	double disij = 0.;
	if (sqrt(rijsq) < 2.05 * a) {
		disij = Distance(p_i, p_j, L, 0);
	} else {
		disij = sqrt(rijsq);
	}
	
	return disij;
}

double Sign(double value, double sign) {
	
	if (fabs(sign) <= 1e-10) {
		return 0.;
	}
	return value * sign / fabs(sign);
}

void IdentifyEnds(Particle *p_i, Particle *p_j, double L, int ends[2]) {
	
	int ellipse1 = 0;
	int ellipse2 = 0;
	
	double rij[dim];
	double rijsq = 0.;
	
	for (int d = 0; d < dim; d++) {
		rij[d] = PeriodicDis(p_j->r[d] - p_i->r[d], L);
		rijsq += rij[d] * rij[d];
	}
	
	double parallel_check = 0;
	
	for (int d = 0; d < dim; d++) {
		
		parallel_check += p_j->e[d] * p_i->e[d];
	}
	
	if (fabs(parallel_check - 1.) < 1e-8) {
		
		//printf("1\n");
		
		double erij[dim];

		for (int d = 0; d < dim; d++) {
			erij[d] = rij[d] / sqrt(rijsq);
		}
		
		double cos_of_angle_ei_erij = 0.;
		for (int d = 0; d < dim; d++) {
			cos_of_angle_ei_erij += p_i->e[d] * erij[d];
		}
		
		if (cos_of_angle_ei_erij > 0) {
			ellipse1 = 1;
			ellipse2 = 0;
		} else {
			ellipse1 = 0;
			ellipse2 = 1;
		}		
	} else if (fabs(- parallel_check - 1.) < 1e-8) {
		
		//printf("2\n");
		
		
		double erij[dim];

		for (int d = 0; d < dim; d++) {
			erij[d] = rij[d] / sqrt(rijsq);
		}
		
		double cos_of_angle_ei_erij = 0.;
		for (int d = 0; d < dim; d++) {
			cos_of_angle_ei_erij += p_i->e[d] * erij[d];
		}
		
		if (cos_of_angle_ei_erij > 0) {
			ellipse1 = 1;
			ellipse2 = 1;
		} else {
			ellipse1 = 0;
			ellipse2 = 0;
		}
	} else {
		
		//printf("3\n");	
		
		double rij_new[dim];
		
		for (int d = 0; d < dim; d++) {
			rij_new[d] = -(a - 0.5) / 2. * p_i->e[d] + rij[d] + (a - 0.5) / 2. * p_j->e[d];
		}
		
		double rijei = 0.;
		double rijej = 0.;
		double eiej = parallel_check;
		double l_half_i = (a - 0.5) / 2.;
		double l_half_j = (a - 0.5) / 2.;
		
		for (int d = 0; d < dim; d++) {
			
			rijei += rij_new[d] * p_i->e[d];
			rijej += rij_new[d] * p_j->e[d];
		}
		
		double denom_temp = 1.0 - eiej * eiej;
		double lambda_prime;
		double mu_prime;

		lambda_prime = (rijei - eiej * rijej) / denom_temp;
		mu_prime = (-rijej + eiej * rijei) / denom_temp;
			
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
		
		//printf("%f\t%f\n",lambda_prime, (a -0.5) / 2.);
		//printf("%f\n",mu_prime);
		
		if (fabs(lambda_prime + (a - 0.5) / 2.) < 1e-8) {
			ellipse1 = 0;
		} else {
			ellipse1 = 1;
		}
		
		if (fabs(mu_prime + (a - 0.5) / 2.) < 1e-8) {
			ellipse2 = 0;
		} else {
			ellipse2 = 1;
		}
	}
	
	ends[0] = ellipse1;
	ends[1] = ellipse2;
}

double Distance(Particle *p_i, Particle *p_j, double L, int state) {
	
	int ends[2];
	IdentifyEnds(p_i, p_j, L, ends);
	
	double a1;
	double b1;	
	if (ends[0] == 1) {
		a1 = p_i->a;
		b1 = p_i->b;
	} else {
		a1 = p_i->b * 1.001;
		b1 = p_i->b;		
	}	
	double ecc1 = sqrt(1. - b1 * b1 / (a1 * a1));
	
	double a2;
	double b2;
	if (ends[1] == 1) {
		a2 = p_j->a;
		b2 = p_j->b;
	} else {
		a2 = p_j->b * 1.001;
		b2 = p_j->b;		
	}
	double ecc2 = sqrt(1. - b2 * b2 / (a2 * a2));
	
	double eta = a1 / b1 - 1.;
	
	double rij[dim];
	double rijsq = 0.;
	
	for (int d = 0; d < dim; d++) {
		rij[d] = PeriodicDis(p_j->r[d] - p_i->r[d], L);
		rijsq += rij[d] * rij[d];
	}
	
	double rijei = 0;
	for (int d = 0; d < dim; d++) {
		
		rijei += rij[d] * p_i->e[d];
	}
	
	double rij_hat[dim];
	double rij_norm = sqrt(rijsq);
	
	for (int d = 0; d < dim; d++) {
		
		rij_hat[d] = rij[d] / rij_norm;
	}
	
	double eiej = 0.;	
	double rij_hatei = 0.;
	double rij_hatej = 0.;
	
	for (int d = 0; d < dim; d++) {
		
		rij_hatei += rij_hat[d] * p_i->e[d];
		rij_hatej += rij_hat[d] * p_j->e[d];
		eiej += p_i->e[d] * p_j->e[d];
	}
	
	double check_para = 1.0 - eiej * eiej;
	if (check_para  < 1e-10) {
		
		if (eiej / fabs(eiej) < 0.) {
			eiej = - eiej;
			rij_hatej = - rij_hatej;
		}
	}	
	
	double d = 0.;
	
	double ratio_bsq = b1 * b1 / (b2 * b2);
	double one_plus_eiej = (1. + eiej);
	double one_minus_eiej = (1. - eiej);
	double ecc1_sq = ecc1 * ecc1;
	double ecc2_sq = ecc2 * ecc2;
	double eiej_sq = eiej * eiej;
	double rij_hat_ei_sq = rij_hatei * rij_hatei;
	double b1_div_a1 = b1 / a1;
	double rij_hat_eiej = rij_hatei * eiej;
	double etaeiej = eta * eiej;
	double eta_two_plus_eta = eta * (2. + eta);
	
	double A11 = ratio_bsq * ( 1. + 0.5 * one_plus_eiej * (eta_two_plus_eta - ecc2_sq * (1. + etaeiej) * (1. + etaeiej)));
	double A22 = ratio_bsq * ( 1. + 0.5 * one_minus_eiej * (eta_two_plus_eta - ecc2_sq * (1. - etaeiej) * (1. - etaeiej)));
	double temp_argument_A12 = check_para;
	if (temp_argument_A12 < 0.) {
		temp_argument_A12 = 0.;
	}	
	double A12 = ratio_bsq * 0.5 * sqrt(temp_argument_A12) * (eta * (2. + eta) + ecc2_sq * (1. - eta * eta * eiej_sq));
	
	double A12_sq = A12 * A12;
	
	double ev_part_1 = 0.5 * (A11 + A22);
	double ev_part_2 = sqrt(0.25 * (A11 - A22) * (A11 - A22) + A12_sq);
	
	double lambda_p = ev_part_1 + ev_part_2;
	double lambda_m = ev_part_1 - ev_part_2;
	
	double e_p_primed_prime = 0.;
	
	if (check_para  < 1e-10) {
		
		if (A11 >= A22) {			
			e_p_primed_prime = 1. / (sqrt(1. - ecc1_sq * rij_hat_ei_sq)) * b1_div_a1 * rij_hatei;
		} else if (A22 < A11) {
			double temp_argument_e_p_primed_prime = 1. - rij_hat_ei_sq;
			if (temp_argument_e_p_primed_prime < 0.) {
				temp_argument_e_p_primed_prime = 0.;
			}
			e_p_primed_prime = sqrt(temp_argument_e_p_primed_prime) / (sqrt(1. - ecc1_sq * rij_hat_ei_sq));
		}
	} else {
		
		double lambda_m_A11 = lambda_p - A11;
		
		double factor = 1. / (sqrt(2.) * sqrt(A12_sq + lambda_m_A11 * lambda_m_A11) * sqrt(1. - ecc1_sq * rij_hat_ei_sq));
		double term1 = A12 / (sqrt(one_plus_eiej)) * (b1_div_a1 * rij_hatei + rij_hatej + (b1_div_a1 - 1.) * rij_hat_eiej);
		double term2 = (lambda_p - A11) / (sqrt(one_minus_eiej)) * (b1_div_a1 * rij_hatei - rij_hatej - (b1_div_a1 - 1.) * rij_hat_eiej);
		
		e_p_primed_prime = factor * (term1 + term2);
	}
	
	double a2_prime = 1. / sqrt(lambda_m);
	double b2_prime = 1. / sqrt(lambda_p);
	
	double delta_dis = a2_prime * a2_prime / (b2_prime * b2_prime) - 1.;
	
	double one_plus_delta = (1. + delta_dis);
	
	double d_prime = 0.;
	
	if (fabs(e_p_primed_prime) < 1e-8 || fabs(delta_dis) < 1e-8) {
		
		d_prime = 1. + a2_prime;
	} else if (fabs(fabs(e_p_primed_prime) - 1.) < 1e-8 && state == 0) { 
		return Distance(p_j, p_i, L, 1);
	} else {
		
		double tanphisq = 1. / (e_p_primed_prime * e_p_primed_prime) - 1.;
		
		double iv_b2_prime = 1. / b2_prime;
		double iv_b2_prime_sq = 1. / (b2_prime * b2_prime);
		double one_plus_tan_plus_delta = (1. + tanphisq + delta_dis);
		double one_plus_tan = (1. + tanphisq);
		
		double A = - iv_b2_prime_sq * one_plus_tan;
		double B = - 2. * iv_b2_prime * one_plus_tan_plus_delta;
		double C = - tanphisq - one_plus_delta * one_plus_delta + iv_b2_prime_sq * (1. + one_plus_delta * tanphisq);
		double D = 2. * iv_b2_prime * one_plus_tan * one_plus_delta;
		double E = one_plus_tan_plus_delta * one_plus_delta;
		
		double A_sq = A * A;
		double B_sq = B * B;
		
		double alpha = - 3. * B_sq/ (8. * A_sq) + C / A;
		double beta = B * B_sq / (8. * A * A_sq) - B * C / (2. * A_sq) + D / A;
		double gamma = -3. * B_sq * B_sq / (256. * A_sq * A_sq) + C * B_sq / (16. * A * A_sq) - B * D / (4. * A_sq) + E / A;
		
		double alpha_sq = alpha * alpha;
		
		double q = 0.;
		if (fabs(beta) > 1e-8) {				
			
			double P = - alpha_sq / 12. - gamma;
			double Q = - alpha * alpha_sq / 108. + alpha * gamma / 3. - beta * beta / 8.;
			double temp_argument_Q = Q * Q / 4. + P * P * P /27.;
			if(fabs(temp_argument_Q) < 1e-8 && temp_argument_Q < 0) {
				temp_argument_Q = - temp_argument_Q;
			}
			double U = cbrt(- Q / 2. + sqrt(temp_argument_Q));
			
			double y = 0.;
			if (fabs(U) > 1e-8) {
				y = - 5. / 6. * alpha + U - P / (3. * U);
			} else {
				y = - 5. / 6. * alpha - cbrt(Q);
			}
			q = - B / (4. * A) + 0.5 * (sqrt(alpha + 2. * y) + sqrt(-(3. * alpha + 2. * y + 2. * beta / (sqrt(alpha + 2. * y)))));
		} else {
			q = - B / (4. * A) + sqrt(0.5 * (- alpha + sqrt(alpha_sq - 4. * gamma)));
		}
		
		double first_pow = 1. + (b2_prime * one_plus_delta) / q;
		double second_pow = 1. + b2_prime / q;
		double qsq_minus_one = (q * q - 1.);
		
		d_prime = sqrt(qsq_minus_one / delta_dis * first_pow * first_pow + (1. - qsq_minus_one / delta_dis) * second_pow * second_pow);
	}
	
	d = d_prime / (sqrt(1. - ecc1_sq * rij_hat_ei_sq)) * b1;
	
	if(isnan(d) || isinf(d)) {
		printf("ERROR: Unexpected result while computing the distance.");
		exit(1);
	}
	
	return rij_norm - d;
}