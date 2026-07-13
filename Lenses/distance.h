#ifndef DISTANCE_H
#define DISTANCE_H

#include "NpT_2D_L.h"

void RotationMatrix(double *old, double angle, double *newV);
double ScalarAngle(double *v1, double *v2);
double PeriodicPos(double dr, double L);
double PeriodicDis(double dr, double L);
double DistanceOverlap(Particle *p_i, Particle *p_j, double L);
void CalculateIntersection(double *center_one, double *center_two, double *inter_one, double *inter_two);
double Distance(Particle *p_i, Particle *p_j, double L);

#endif
