#ifndef PRINTDATA_H
#define PRINTDATA_H

#include "NpT_2D_P.h"

void PrintVMD(FILE *f, Particle *P, int step);
void PrintOvito(FILE *f, Particle *P, int step, double L);

#endif