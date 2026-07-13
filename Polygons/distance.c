#include "distance.h"

// ----------------------------------------------------------------------------------------

// Rotates vectors
void RotationMatrix(double *old, double angle, double *newV) {
	
	double *old_temp = old;
	
	newV[0] = old_temp[0] * cos(angle) - old_temp[1] * sin(angle);
	newV[1] = old_temp[0] * sin(angle) + old_temp[1] * cos(angle);
}

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
	double limit = 1.05 / sin(M_PI / (double)m) * a;
	if(lambda > 1.) {
		limit *= lambda;
	}
	if (sqrt(rijsq) < limit) {
		disij = Distance(p_i, p_j, L);
	} else {
		disij = sqrt(rijsq);
	}
	
	return disij;
}

// ----------------------------------------------------------------------------------------

// Calculates the distance between a point and a line segment.
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

// ----------------------------------------------------------------------------------------

// Function which calculates the projection of the particles on the candidate for the 
// separating axis
void ProjectSAT(double ri[][dim], double *axis, double *min, double *max) {
	
	*min = 0.;
	for(int d = 0; d < dim; d++) {
		*min += ri[0][d] * axis[d];
	}
	*max = *min;
	
	
	for(int k = 1; k < m; k++) {
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

// ----------------------------------------------------------------------------------------

// Implementation of the separating axis theorem
int CheckOverlapSAT(double ri[][dim], double rj[][dim], double directionsi[][dim], double directionj[][dim]) {
	
	double mini = 0.;
	double maxi = 0.;
	double minj = 0.;
	double maxj = 0.;
	for(int k = 0; k < m; k++) {
		double axis[2] = {-directionsi[k][1], directionsi[k][0]};
		ProjectSAT(ri, axis, &mini, &maxi);
		ProjectSAT(rj, axis, &minj, &maxj);
		if ((maxi < minj) || (maxj < mini)) {
			return 0;
		}
	}
	
	for(int k = 0; k < m; k++) {
		double axis[2] = {-directionj[k][1], directionj[k][0]};
		ProjectSAT(ri, axis, &mini, &maxi);
		ProjectSAT(rj, axis, &minj, &maxj);
		if ((maxi < minj) || (maxj < mini)) {
			return 0;
		}
	}
	
	return 1;
}

// ----------------------------------------------------------------------------------------

// Deforms the polygon shapes by scaling on direction byb lambda
void Compress(double ri[][dim], double *e) {

	double ri_temp[m][dim];
	double e_perp[2] = {-e[1], e[0]};
	
	for(int k = 0; k < m; k++) {
		
		double ri_m_e = 0.;
		double ri_m_e_perp = 0.;
		for (int d = 0; d < dim; d++) {
			ri_m_e += ri[k][d] * e[d];
			ri_m_e_perp += ri[k][d] * e_perp[d];
		}
		
		for (int d = 0; d < dim; d++) {
			ri_temp[k][d] = ri_m_e * e[d] * lambda + ri_m_e_perp * e_perp[d];
		}
	}
	
	for(int k = 0; k < m; k++) {
		for(int d = 0; d < dim; d++) {
			ri[k][d] = ri_temp[k][d];
		}
	}
}

// ----------------------------------------------------------------------------------------

// Estimates a distance between the two polygons. This function is not important for 
// calculating the overlap
double Distanceij(double ri[][dim], double cornerj[][dim], double directionj[][dim], double *length) {
	
	double r_temp[dim];
	double corner_temp[dim];
	double direction_temp[dim];
	double length_temp = 0.;
	double dis_temp = 0.;
	double dij_temp[dim];
	double dis = 0.;
	
	for (int k = 0; k < m; k++) {
		for (int h = 0; h < m; h++) {
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

// ----------------------------------------------------------------------------------------

// Checks for overlap between two polygons
double Distance(Particle *p_i, Particle *p_j, double L) {
	
	double rij[dim];
	for (int d = 0; d < dim; d++) {
		rij[d] = PeriodicDis(p_j->r[d] - p_i->r[d], L);
	}
	
	double ri[m][dim];	
	double cornerj[m][dim];
	double directioni[m][dim];
	double directionj[m][dim];
	double length = 0.5 / sin(M_PI / (double)m) * a;
	double length_k[m];

	// Constructs the regular polygons via the vectors pointing from the center 
	// to the vertices
	double vector_to_corner_i[m][dim];
	double vector_to_corner_j[m][dim];
	
	for (int d = 0; d < dim; d++) {
		vector_to_corner_i[0][d] = length * p_i->e[d];
		vector_to_corner_j[0][d] = length * p_j->e[d];
	}
	
	for (int i = 1; i < m; i++) {
		RotationMatrix(vector_to_corner_i[i - 1], 2. * M_PI / (double)m, vector_to_corner_i[i]);
		RotationMatrix(vector_to_corner_j[i - 1], 2. * M_PI / (double)m, vector_to_corner_j[i]);
	}

	// Deforms the particles
	if(m % 2 != 0) {	
		Compress(vector_to_corner_i, p_i->e);
		Compress(vector_to_corner_j, p_j->e);
	} else {
		double e_perp_i[2] = {-p_i->e[1], p_i->e[0]};
		double e_perp_j[2] = {-p_j->e[1], p_j->e[0]};
		Compress(vector_to_corner_i, e_perp_i);
		Compress(vector_to_corner_j, e_perp_j);
	}

	// Vectors from the system origin to the vertices
	for (int k = 0; k < m; k++) {                  
		for (int d = 0; d < dim; d++) {
			ri[k][d] = p_i->r[d] + vector_to_corner_i[k][d];
			cornerj[k][d] = p_i->r[d] + rij[d] + vector_to_corner_j[k][d];
		}
	}	

	// Calculates vectors parallel to the polygon edges
	for(int k = 0; k < m; k++) {
		int kp1 = k + 1;
		if (kp1 == m) {
			kp1 = 0;
		}
		double length_k_temp = 0.;
		for (int d = 0; d < dim; d++) {
			directioni[k][d] = ri[kp1][d] - ri[k][d];
			directionj[k][d] = cornerj[kp1][d] - cornerj[k][d];
			length_k_temp += directioni[k][d] * directioni[k][d];
		}
		for (int d = 0; d < dim; d++) {
			directioni[k][d] /= sqrt(length_k_temp);
			directionj[k][d] /= sqrt(length_k_temp);
			length_k[k] = sqrt(length_k_temp);
		}
	}

	// Uses separating axis theorm check for overlap
	int checkOverlap = CheckOverlapSAT(ri, cornerj, directioni, directionj);
	if (checkOverlap == 1) {		
		return -1.;
	}

	// Estimates approximations for the distance if the particle do not 
	// overlap
	double dis1 = Distanceij(ri, cornerj, directionj, length_k);
	double dis2 = Distanceij(cornerj, ri, directioni, length_k);
	
	if (dis1 < dis2) {
		return dis1;
	} else {
		return dis2;
	}
}

// ----------------------------------------------------------------------------------------
