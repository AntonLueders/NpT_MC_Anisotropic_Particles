#ifndef DISTANCE_H
#define DISTANCE_H

#include "NpT_2D_P.h"

double PeriodicPos(double dr, double L);
double PeriodicDis(double dr, double L);
double DistanceOverlap(Particle *p_i, Particle *p_j, double L);
double Sign(double value, double sign);
void IdentifyEnds(Particle *p_i, Particle *p_j, double L, int ends[2]);
double Distance(Particle *p_i, Particle *p_j, double L, int state);

#endif
