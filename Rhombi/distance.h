#ifndef DISTANCE_H
#define DISTANCE_H

#include "NpT_2D_D.h"

void RotationMatrix(double *old, double angle, double *newV);
double PeriodicPos(double dr, double L);
double PeriodicDis(double dr, double L);
double DistanceOverlap(Particle *p_i, Particle *p_j, double L);
double DistancePointLine(double *r, double *corner, double *direction, double length, double *dij);
void ProjectSAT(double ri[][dim], double *axis, double *min, double *max);
int CheckOverlapSAT(double ri[][dim], double rj[][dim], double directionsi[][dim], double directionj[][dim]);
double Distanceij(double ri[][dim], double cornerj[][dim], double directionj[][dim], double *length);
double Distance(Particle *p_i, Particle *p_j, double L);


#endif
