#include "distance.h"

// ----------------------------------------------------------------------------------------

void RotationMatrix(double *old, double angle, double *newV) {
	
	double *old_temp = old;
	
	newV[0] = old_temp[0] * cos(angle) - old_temp[1] * sin(angle);
	newV[1] = old_temp[0] * sin(angle) + old_temp[1] * cos(angle);
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
	if (sqrt(rijsq) < 1.05 * a) {
		disij = Distance(p_i, p_j, L);
	} else {
		disij = sqrt(rijsq);
	}
	
	return disij;
}


// ----------------------------------------------------------------------------------------

double DistancePointLine(double *r, double *corner, double *direction, double length, double *dij) {
	
	double position[dim];
	double parameter = 0.;
	
	for (int d = 0; d < dim; d++) {
		position[d] = r[d] - corner[d];
		parameter += position[d] * direction[d] / length;
	}
	
	double point_on_line[dim];
	
	if (parameter < 0.) {
		for (int d = 0; d < dim; d++) {
			point_on_line[d] = corner[d];
		}
	} else if (parameter > 1.) {
		for (int d = 0; d < dim; d++) {
			point_on_line[d] = corner[d] + length * direction[d];
		}
	} else {
		for (int d = 0; d < dim; d++) {
			point_on_line[d] = corner[d] + parameter * length * direction[d];
		}
	}
	
	double dijsq = 0.;
	
	for (int d = 0; d < dim; d++) {
		dij[d] = point_on_line[d] - r[d];
		dijsq += dij[d] * dij[d];
	}
	
	return sqrt(dijsq);
}

void ProjectSAT(double ri[][dim], double *axis, double *min, double *max) {
	
	*min = 0.;
	for(int d = 0; d < dim; d++) {
		*min += ri[0][d] * axis[d];
	}
	*max = *min;
	
	
	for(int k = 1; k < 4; k++) {
		double scalar = 0.;
		for(int d = 0; d < dim; d++) {
			scalar += ri[k][d] * axis[d];
		}
		if (scalar < *min) {
			*min = scalar;
		}
		if (scalar > *max) {
			*max = scalar;
		}
	}
}

int CheckOverlapSAT(double ri[][dim], double rj[][dim], double directionsi[][dim], double directionj[][dim]) {
	
	double mini = 0.;
	double maxi = 0.;
	double minj = 0.;
	double maxj = 0.;
	for(int k = 0; k < 4; k++) {
		double axis[2] = {-directionsi[k][1], directionsi[k][0]};
		ProjectSAT(ri, axis, &mini, &maxi);
		ProjectSAT(rj, axis, &minj, &maxj);
		if ((maxi < minj) || (maxj < mini)) {
			return 0;
		}
	}
	
	for(int k = 0; k < 4; k++) {
		double axis[2] = {-directionj[k][1], directionj[k][0]};
		ProjectSAT(ri, axis, &mini, &maxi);
		ProjectSAT(rj, axis, &minj, &maxj);
		if ((maxi < minj) || (maxj < mini)) {
			return 0;
		}
	}
	
	return 1;
}

double Distanceij(double ri[][dim], double cornerj[][dim], double directionj[][dim], double *length) {
	
	double r_temp[dim];
	double corner_temp[dim];
	double direction_temp[dim];
	double length_temp = 0.;
	double dis_temp = 0.;
	double dij_temp[dim];
	double dis = 0.;
	
	for (int k = 0; k < 4; k++) {
		for (int h = 0; h < 4; h++) {
			for(int d = 0; d < dim; d++) {
				r_temp[d] = ri[k][d];
				corner_temp[d] = cornerj[h][d];
				direction_temp[d] = directionj[h][d];
				length_temp = length[h];
			}
			
			dis_temp = DistancePointLine(r_temp, corner_temp, direction_temp, length_temp, dij_temp);
			if (k == 0 && h == 0) {
				dis = fabs(dis_temp);
			}
			
			if (dis_temp < dis) {
				dis = fabs(dis_temp);
			}
		}		
	}
	
	return dis;
}

double Distance(Particle *p_i, Particle *p_j, double L) {
	
	double rij[dim];
	for (int d = 0; d < dim; d++) {
		rij[d] = PeriodicDis(p_j->r[d] - p_i->r[d], L);
	}
	
	double ri[4][dim];	
	double cornerj[4][dim];
	double directioni[4][dim];
	double directionj[4][dim];
	double length[4];
	double length_k[4];
	
	for(int k = 0; k < 4; k++) {
		if (k % 2 == 0) {
			length[k] = a / 2.;
		} else {
			length[k] = b / 2.;
		}
	}
	
	double vector_to_corner_i[4][dim];
	double vector_to_corner_j[4][dim];
	
	for (int d = 0; d < dim; d++) {
		vector_to_corner_i[0][d] = p_i->e[d];
		vector_to_corner_j[0][d] = p_j->e[d];
	}
	
	for (int i = 1; i < 4; i++) {
		RotationMatrix(vector_to_corner_i[i - 1], 2. * M_PI / 4., vector_to_corner_i[i]);
		RotationMatrix(vector_to_corner_j[i - 1], 2. * M_PI / 4., vector_to_corner_j[i]);
	}
	
	for (int k = 0; k < 4; k++) {                  
		for (int d = 0; d < dim; d++) {
			ri[k][d] = p_i->r[d] + vector_to_corner_i[k][d] * length[k];
			cornerj[k][d] = p_i->r[d] + rij[d] + vector_to_corner_j[k][d] * length[k];
		}
		length_k[k] = sqrt((a / 2.) * (a / 2.) + (b / 2.) * (b / 2.));
	}
		
	for(int k = 0; k < 4; k++) {
		int kp1 = k + 1;
		if (kp1 == 4) {
			kp1 = 0;
		}
		for (int d = 0; d < dim; d++) {
			directioni[k][d] = (ri[kp1][d] - ri[k][d]) / length_k[k];
			directionj[k][d] = (cornerj[kp1][d] - cornerj[k][d]) / length_k[k];
		}
	}
	
	int checkOverlap = CheckOverlapSAT(ri, cornerj, directioni, directionj);
	if (checkOverlap == 1) {
		return -1.;
	}
	
	double dis1 = Distanceij(ri, cornerj, directionj, length_k);
	double dis2 = Distanceij(cornerj, ri, directioni, length_k);
	
	if (dis1 < dis2) {
		return dis1;
	} else {
		return dis2;
	}
}
