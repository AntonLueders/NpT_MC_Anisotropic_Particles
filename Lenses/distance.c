#include "distance.h"

// ----------------------------------------------------------------------------------------

void RotationMatrix(double *old, double angle, double *newV) {
	
	double *old_temp = old;
	
	newV[0] = old_temp[0] * cos(angle) - old_temp[1] * sin(angle);
	newV[1] = old_temp[0] * sin(angle) + old_temp[1] * cos(angle);
}

double ScalarAngle(double *v1, double *v2) {
	
	double result = 0.;
	double abs_v1 = 0.; 
	double abs_v2 = 0.;
	
	for (int d = 0; d < dim; d++) {
		result += v1[d] * v2[d];
		abs_v1 += v1[d] * v1[d];
		abs_v2 += v2[d] * v2[d];
	}
	
	return acos(result / sqrt(abs_v1) / sqrt(abs_v2));
}

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
	double limit = 2.05 * a;
	if (sqrt(rijsq) < limit) {
		disij = Distance(p_i, p_j, L);
	} else {
		disij = sqrt(rijsq);
	}
	
	return disij;
}

void CalculateIntersection(double *center_one, double *center_two, double *inter_one, double *inter_two) {
	
	double diss = 0.;
	double diss_vec[dim];
	
	for (int d = 0; d < dim; d++) {
		diss_vec[d] = center_two[d] - center_one[d];
		diss += diss_vec[d] * diss_vec[d];
	}
	
	diss = sqrt(diss);
	
	if (diss < pow(10., -8.)) {
		printf("ERROR: Circle intersection error.");
        exit(1);
	}
	
	double x = diss / 2.;
	
	double arg_y = R * R - x * x;
	if (arg_y < pow(10., -8.)) {
		arg_y = 0.;
	}
	
	double y = sqrt(arg_y);
	
	inter_one[0] = center_one[0] + x * diss_vec[0] / diss - y * diss_vec[1] / diss;
	inter_one[1] = center_one[1] + x * diss_vec[1] / diss + y * diss_vec[0] / diss;
	
	inter_two[0] = center_one[0] + x * diss_vec[0] / diss + y * diss_vec[1] / diss;
	inter_two[1] = center_one[1] + x * diss_vec[1] / diss - y * diss_vec[0] / diss;
}

double Distance(Particle *p_i, Particle *p_j, double L) {
	
	double rij[dim];
	double diss_rij = 0.;
	for (int d = 0; d < dim; d++) {
		rij[d] = PeriodicDis(p_j->r[d] - p_i->r[d], L);
		diss_rij += rij[d] * rij[d];
	}
	
	diss_rij = sqrt(diss_rij);
	
	if(diss_rij <= 2. * b) {
		return -1.;
	}
	
	double alpha = asin(a / R);
	
	double e_perp_i[dim];
	double e_perp_j[dim];
	
	e_perp_i[0] = -p_i->e[1];
	e_perp_i[1] = p_i->e[0];
	
	e_perp_j[0] = -p_j->e[1];
	e_perp_j[1] = p_j->e[0];
	
	double c_i[2][dim];
	double c_j[2][dim];
	
	for (int d = 0; d < dim; d++) {
		c_i[0][d] = p_i->r[d] - (R - b) * e_perp_i[d];
		c_i[1][d] = p_i->r[d] + (R - b) * e_perp_i[d];
		c_j[0][d] = p_i->r[d] + rij[d] - (R - b) * e_perp_j[d];
		c_j[1][d] = p_i->r[d] + rij[d] + (R - b) * e_perp_j[d];
	}
	
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			
			double c1[dim];
			double c2[dim];
			
			for (int d = 0; d < dim; d++) {
				
				c1[d] = c_i[i][d];
				c2[d] = c_j[j][d];
			}
			
			double inter_one[dim];
			double inter_two[dim];
			
			double diss_c1_c2 = 0.;
			for (int d = 0; d < dim; d++) {
				diss_c1_c2 += (c1[d] - c2[d]) * (c1[d] - c2[d]);
			}
			
			diss_c1_c2 = sqrt(diss_c1_c2);
			
			if(diss_c1_c2 > 2. * R) {
				continue;
			}
			
			CalculateIntersection(c1, c2, inter_one, inter_two);
			
			double Delta_c1[dim];
			double Delta_inter_one_c1[dim];
			double Delta_inter_two_c1[dim];
			double Delta_c2[dim];
			double Delta_inter_one_c2[dim];
			double Delta_inter_two_c2[dim];
			
			for (int d = 0; d < dim; d++) {
				Delta_c1[d] = p_i->r[d] - c1[d];
				Delta_c2[d] = p_i->r[d] + rij[d] - c2[d];
				
				Delta_inter_one_c1[d] = inter_one[d] - c1[d];
				Delta_inter_two_c1[d] = inter_two[d] - c1[d];
				
				Delta_inter_one_c2[d] = inter_one[d] - c2[d];
				Delta_inter_two_c2[d] = inter_two[d] - c2[d];
			}
			
			double angle_11 = ScalarAngle(Delta_c1, Delta_inter_one_c1);
			double angle_12 = ScalarAngle(Delta_c1, Delta_inter_two_c1);
			double angle_21 = ScalarAngle(Delta_c2, Delta_inter_one_c2);
			double angle_22 = ScalarAngle(Delta_c2, Delta_inter_two_c2);
			
			if ((angle_11 <= alpha && angle_21 <= alpha) || (angle_12 <= alpha && angle_22 <= alpha)) {
				return -1.;
			}
		}
	}
	
	return 2. * b;
}