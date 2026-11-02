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
	
	if(sqrt(rijsq) < 2. * b) {
		return sqrt(rijsq) - 2. * b;
	}
	
	double disij = 0.;
	if (sqrt(rijsq) < 2.05 * a) {
		disij = Distance(p_i, p_j, L, 0);
	} else {
		disij = sqrt(rijsq) - 2. * a;
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

// This function identifies which ends are interacting to properly assign either the disk
// or the ellipse shape
int IdentifyEnds(Particle *p_i, Particle *p_j, double L, int ends[2]) {
	
	int ellipse1 = 0;
	int ellipse2 = 0;
	
	double rij[dim];
	double rijsq = 0.;
	
	for (int d = 0; d < dim; d++) {
		rij[d] = PeriodicDis(p_j->r[d] - p_i->r[d], L);
		rijsq += rij[d] * rij[d];
	}
	
	//Check whether the two internal egg line segments intersect.
	{
		double segment_length = a - b;

		// Start point of segment i is chosen as the local origin.
		double qx = rij[0];
		double qy = rij[1];

		// Direction vectors of the two finite internal segments.
		double rx = segment_length * p_i->e[0];
		double ry = segment_length * p_i->e[1];

		double sx = segment_length * p_j->e[0];
		double sy = segment_length * p_j->e[1];

		double cross_rs = rx * sy - ry * sx;
		double cross_qr = qx * ry - qy * rx;

		if (fabs(cross_rs) > 1e-12) {

			// Segment intersection parameters:
			// p_i->r + t r = p_i->r + rij + u s
			double t = (qx * sy - qy * sx) / cross_rs;
			double u = cross_qr / cross_rs;

			if (t >= -1e-12 && t <= 1. + 1e-12 &&
			    u >= -1e-12 && u <= 1. + 1e-12) {
				return 0;
			}

		} else {

			// Parallel segments.
			double cross_qs = qx * sy - qy * sx;

			if (fabs(cross_qs) <= 1e-12) {

				double r_sq = rx * rx + ry * ry;

				double t0 = (qx * rx + qy * ry) / r_sq;
				double t1 = ((qx + sx) * rx + (qy + sy) * ry) / r_sq;

				double t_min = fmin(t0, t1);
				double t_max = fmax(t0, t1);

				if (t_max >= -1e-12 && t_min <= 1.0 + 1e-12) {
					return 0;
				}
			}
		}
	}
	
	double parallel_check = 0;
	
	for (int d = 0; d < dim; d++) {
		
		parallel_check += p_j->e[d] * p_i->e[d];
	}
	
	if (fabs(parallel_check - 1.) < 1e-8) {

		double erij[dim];

		for (int d = 0; d < dim; d++) {
			erij[d] = rij[d] / sqrt(rijsq);
		}
		
		double cos_of_angle_ei_erij = 0.;
		for (int d = 0; d < dim; d++) {
			cos_of_angle_ei_erij += p_i->e[d] * erij[d];
		}
		
		if (cos_of_angle_ei_erij >= 0) {
			ellipse1 = 1;
			ellipse2 = 0;
		} else {
			ellipse1 = 0;
			ellipse2 = 1;
		}		
	} else if (fabs(- parallel_check - 1.) < 1e-8) {
		
		double erij[dim];

		for (int d = 0; d < dim; d++) {
			erij[d] = rij[d] / sqrt(rijsq);
		}
		
		double cos_of_angle_ei_erij = 0.;
		for (int d = 0; d < dim; d++) {
			cos_of_angle_ei_erij += p_i->e[d] * erij[d];
		}
		
		if (cos_of_angle_ei_erij >= 0) {
			ellipse1 = 1;
			ellipse2 = 1;
		} else {
			ellipse1 = 0;
			ellipse2 = 0;
		}
	} else {
		
		double rij_new[dim];
		
		for (int d = 0; d < dim; d++) {
			rij_new[d] = -(a - b) / 2. * p_i->e[d] + rij[d] + (a - b) / 2. * p_j->e[d];
		}
		
		double rijei = 0.;
		double rijej = 0.;
		double eiej = parallel_check;
		double l_half_i = (a - b) / 2.;
		double l_half_j = (a - b) / 2.;
		
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
		
		if (fabs(lambda_prime + (a - b) / 2.) < 1e-8) {
			ellipse1 = 0;
		} else {
			ellipse1 = 1;
		}
		
		if (fabs(mu_prime + (a - b) / 2.) < 1e-8) {
			ellipse2 = 0;
		} else {
			ellipse2 = 1;
		}
	}
	
	ends[0] = ellipse1;
	ends[1] = ellipse2;
	
	return 1;
}

// ----------------------------------------------------------------------------------------

// Contact distance between two ellipses
// Implemented following the supporting informations of
// Zheng, X. & Palffy-Muhoray, P. Distance of closest approach of two arbitrary hard ellipses in two dimensions.
// Phys. Rev. E 75, 061709; 10.1103/PhysRevE.75.061709 (2007).
double DistanceEllipses(Particle *p_i, Particle *p_j, double L, int state) {
	
	double a1 = p_i->a;
	double b1 = p_i->b;
	double ecc1 = p_i->ecc;
	
	double b2 = p_j->b;
	double ecc2 = p_j->ecc;
	
	double eta = a1 / b1 - 1.;
	
	double rij[dim];
	double rijsq = 0.;
	
	for (int d = 0; d < dim; d++) {
		rij[d] = PeriodicDis(p_j->r[d] - p_i->r[d], L);
		rijsq += rij[d] * rij[d];
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
	if (eiej > 1.0) {
		eiej = 1.0;
	} else if (eiej < -1.0) {
		eiej = -1.0;
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
		} else if (A22 > A11) {
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
			double temp_argument_Q = Q * Q / 4. + P * P * P / 27.;
			if (temp_argument_Q < 0.0) {
				temp_argument_Q = 0.;
			}
			double U = cbrt(- Q / 2. + sqrt(temp_argument_Q));
			
			double y = 0.;
			if (fabs(U) > 1e-8) {
				y = - 5. / 6. * alpha + U - P / (3. * U);
			} else {
				y = - 5. / 6. * alpha - cbrt(Q);
			}
			
			double alpha_plus_2y = alpha + 2. * y;			
			if (alpha_plus_2y < 0.) {
				alpha_plus_2y = 0.;
			}

			double minus_3_alpha_plus_2y_plus_2beta_norm = -(3. * alpha + 2. * y + 2. * beta / (sqrt(alpha_plus_2y) + 1e-14));
			if (minus_3_alpha_plus_2y_plus_2beta_norm < 0.) {
				minus_3_alpha_plus_2y_plus_2beta_norm = 0.;
			}
			
			q = - B / (4. * A) + 0.5 * (sqrt(alpha_plus_2y) + sqrt(minus_3_alpha_plus_2y_plus_2beta_norm));
		} else {

			double alpha_sq_minus_4_gamma = alpha_sq - 4. * gamma;
			if (alpha_sq_minus_4_gamma < 0.) {
				alpha_sq_minus_4_gamma = 0.;
			}

			double half_of_minus_alpha_plus_root = 0.5 * (- alpha + sqrt(alpha_sq_minus_4_gamma));
			if (half_of_minus_alpha_plus_root < 0.) {
				half_of_minus_alpha_plus_root = 0.;
			}
			
			q = - B / (4. * A) + sqrt(half_of_minus_alpha_plus_root);
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

// ----------------------------------------------------------------------------------------

// Distance between disk-shaped ends
double DistanceDisks(Particle *p_i, Particle *p_j, double L) {
	
	double rijsq = 0.;
	for (int d = 0; d < dim; d++) {
		rijsq += PeriodicDis(p_j->r[d] - p_i->r[d], L) * PeriodicDis(p_j->r[d] - p_i->r[d], L);
	}
	
	return sqrt(rijsq) - 2. * b;
}

// ----------------------------------------------------------------------------------------

// Signed surface distance between a disk and an ellipse.
double DistanceDiskEllipse(Particle *p_circle, Particle *p_ellipse, double L) {

    // Tolerances
    double tol_F = 1e-8;         	// residual tolerance on F(lambda) 
    double tol_lambda_rel = 1e-8;	// relative tolerance on bracket width 
    int    max_iter = 30;         	// iterations of hybrid solver 

    double R = p_circle->b;

    double a_ellipse = p_ellipse->a;
    double b_ellipse = p_ellipse->b;

    double rx = PeriodicDis(p_circle->r[0] - p_ellipse->r[0], L);
    double ry = PeriodicDis(p_circle->r[1] - p_ellipse->r[1], L);

    // Ellipse body-frame basis. 
    double ex = p_ellipse->e[0];
    double ey = p_ellipse->e[1];
    double epx = -ey;
    double epy = ex;

    // Disk-center coordinates in ellipse frame. 
    double x = rx * ex + ry * ey;
    double y = rx * epx + ry * epy;

    double a_sq = a_ellipse * a_ellipse;
    double b_sq = b_ellipse * b_ellipse;

    double ellipse_coordinate = x * x / a_sq + y * y / b_sq;
    if (ellipse_coordinate <= 1.) {
        return -R;
    }
	
    double lambda_lo = 0.;

    // Find an upper bracket lambda_hi with F(lambda_hi) <= 0.
    double lambda_hi = a_sq + b_sq; 		// reasonable starting scale 
    double denom_a_hi = a_sq + lambda_hi;
    double denom_b_hi = b_sq + lambda_hi;

    double F_hi = a_sq * x * x / (denom_a_hi * denom_a_hi)
                + b_sq * y * y / (denom_b_hi * denom_b_hi)
                - 1.;

    int guard = 0;
    while (F_hi > 0.) {
        lambda_hi *= 2.;
        denom_a_hi = a_sq + lambda_hi;
        denom_b_hi = b_sq + lambda_hi;

        F_hi = a_sq * x * x / (denom_a_hi * denom_a_hi)
             + b_sq * y * y / (denom_b_hi * denom_b_hi)
             - 1.;

        if (++guard > 100) {
            printf("ERROR: Could not bracket disk-ellipse closest-point root.\n");
            exit(EXIT_FAILURE);
        }
    }

    // Start hybrid solve at the mid-bracket
    double lambda = 0.5 * (lambda_lo + lambda_hi);

    for (int iter = 0; iter < max_iter; iter++) {

        double denom_a = a_sq + lambda;
        double denom_b = b_sq + lambda;

        double Fa = a_sq * x * x / (denom_a * denom_a);
        double Fb = b_sq * y * y / (denom_b * denom_b);

        double F = Fa + Fb - 1.;

        if (fabs(F) < tol_F) {
            break;
        }
        if ((lambda_hi - lambda_lo) <= tol_lambda_rel * (1.0 + lambda_hi)) {
            break;
        }

        // Derivative: F'(lambda) < 0 for lambda >= 0
        double Fp = -2.0 * a_sq * x * x / (denom_a * denom_a * denom_a)
                    -2.0 * b_sq * y * y / (denom_b * denom_b * denom_b);

        // Update bracket using F sign at current lambda 
        if (F > 0.) {
            // lambda too small 
            lambda_lo = lambda;
        } else {
            // lambda too large 
            lambda_hi = lambda;
        }

        // Attempt a Newton step, safeguarded by the bracket 
        double take_newton = 1.;
        double lambda_new = lambda;

        if (!(isfinite(Fp)) || Fp >= 0.) {
            take_newton = 0;
        } else {
            double step = F / Fp;  // Fp < 0, so this moves the right way 
            lambda_new = lambda - step;

            if (!(isfinite(lambda_new)) ||
                lambda_new <= lambda_lo ||
                lambda_new >= lambda_hi) {
                take_newton = 0;
            }
        }

        if (!take_newton) {
            // Fall back to bisection for this iteration
            lambda_new = 0.5 * (lambda_lo + lambda_hi);
        }

        lambda = lambda_new;
    }

    // Closest point on ellipse boundary (ellipse frame) 
    double X = a_sq * x / (a_sq + lambda);
    double Y = b_sq * y / (b_sq + lambda);

    // Distance from disk center to ellipse boundary
    double dx = x - X;
    double dy = y - Y;
    double point_to_ellipse_distance = sqrt(dx * dx + dy * dy);

    // Contact distance
    return point_to_ellipse_distance - R;
}

// ----------------------------------------------------------------------------------------

double Distance(Particle *p_i, Particle *p_j, double L, int state) {
	
	int ends[2];
	int check_no_inter = IdentifyEnds(p_i, p_j, L, ends);

	if (!check_no_inter) {
		return -1.;
	}
	
	if (ends[0] == 0 && ends[1] == 0) {
		//printf("CircleCircle\n");
		return DistanceDisks(p_i, p_j, L);
	}

	if (ends[0] == 0 && ends[1] == 1) {
		//printf("CircleEps\n");
		return DistanceDiskEllipse(p_i, p_j, L);
	}

	if (ends[0] == 1 && ends[1] == 0) {
		//printf("EpsCircle\n");
		return DistanceDiskEllipse(p_j, p_i, L);
	}

	if (ends[0] == 1 && ends[1] == 1) {
		//printf("EpsEps\n");
		return DistanceEllipses(p_i, p_j, L, state);
	}
	
	return -1.;
}

// ----------------------------------------------------------------------------------------
