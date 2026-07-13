#ifndef NpT_2D_9_H
#define NpT_2D_9_H

#include "lists.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_log.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define DIM 2
#define NP 500
#define PRESSURE 3.
#define TEMP 1.
#define SWEEPS 100000
#define RATE 100
#define AA 0.75
#define BB 0.5
#define SEED 1
#define DELTA 0.17
#define DELTAR 0.53
#define DELTAV 1.1
#define RV 1.55
#define M 3
#define LAMBDA 1.0

// Definition of bool-variables
typedef enum { false, true } bool;

// Structure which manages the particle information
struct particle {
    double r[DIM]; // Position 
    double r_old[DIM];
	double e[DIM];
	double e_old[DIM];
    double a;
	double b;
    Node *verList; // Verlet-list
};
typedef struct particle Particle;

extern gsl_rng *generator;
extern int dim;
extern int N;
extern double p;
extern double T;
extern int sweeps;
extern int rate;
extern double a;
extern double b;
extern long seed;
extern double delta;
extern double delta_r;
extern double delta_V;
extern double r_v;
extern double R;

#endif