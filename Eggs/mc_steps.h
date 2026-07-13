#ifndef MC_STEPS_H
#define MC_STEPS_H

#include "NpT_2D_P.h"

int MetropolisTranslation(int index, Particle *P, double L);
int MetropolisRotation(int index, Particle *P, double L);
double VolumeMCStep(double L, Particle *P);

#endif