#ifndef DISTANCE_H
#define DISTANCE_H

#include "NpT_2D_E.h"

double PeriodicPos(double dr, double L);
double PeriodicDis(double dr, double L);
double DistanceOverlap(Particle *p_i, Particle *p_j, double L);
double Distance(Particle *p_i, Particle *p_j, double L, int state);

#endif
