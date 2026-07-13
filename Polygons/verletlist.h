#ifndef VERLETLIST_H
#define VERLETLIST_H

#include "NpT_2D_M.h"
#include "lists.h"

void GenerateVerletLists(Particle *P, double L, int step);
bool CheckVerlet(Particle *P, double L);
void FreeVerLists(Particle *P);

#endif
