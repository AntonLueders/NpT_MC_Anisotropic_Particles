#include "NpT_2D_M.h"
#include "init.h"
#include "mc_steps.h"
#include "printdata.h"
#include "distance.h"
#include "verletlist.h"
#include "inputdata.h"
#include "printdisplay.h"

// ----------------------------------------------------------------------------------------
// Global variables
// ----------------------------------------------------------------------------------------    

gsl_rng *generator = NULL;                                  // Random number generator
int dim = DIM;
int N = NP;
double p = PRESSURE;
double T = TEMP;
int sweeps = SWEEPS;
int rate = RATE;
double a = AA;
long seed = SEED;
double delta = DELTA;
double delta_r = DELTAR;
double delta_V = DELTAV;
double r_v = RV;
int m = M;
double lambda = LAMBDA;

// ----------------------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------------------- 

int main(int argc, char *argv[]) {
	
	// For the calculation of the total simulation time. 
    time_t start;
    time(&start);
	
	PrintWelcome();
	
	// The program only runs with "inputfile" in the same folder. 
    // A new inputfile can be generated with "./NpT_2D_M.out Setup"
    if (argc == 2 && !strcmp(argv[1], "Setup")) {

        // Generates an example input file which can be customized. 
        PrintInputfile();     // See inputdata.c

        time_t now;
        time(&now);

        // Total simulation time.
        long simTime = now - start;
        printf("Duration of the simulation: %ld s\n", simTime);

        return 0;
    }

	// Setup 
	ReadInput();
	PrintInput();
	CheckForConflict();

	// Sets initial system area
	double V = (double)N * T / p;
	double V_init_min = 6.0 * (double)N * T / 3.;
	if (V < V_init_min) {
		V = V_init_min;
	}
	double L = pow(V, 1. / (double)dim);
	
	printf("init. V: %f\t init. L: %f\n", V, L);
	
	Particle P[N];
	
	// See GSL documentation for more info regarding the generator.
    // GLS has to be linked to get random numbers!
    generator = gsl_rng_alloc(gsl_rng_mt19937);
    printf("Random number generator: mt19937\n");
    gsl_rng_set(generator, seed);

	// Initializes particle positions and particle orientations randomly
	InitRandom(P, L);
	
	double acc_rate[10][3];
	long counter_verlet = 1;
	long counter_acc_translation = 0;
	long counter_acc_rotation = 0;
	long counter_acc_volume = 0;
	long counter_store_current_rate = 0;
	
	FILE *position_file_ovito = NULL;
    position_file_ovito = fopen("Position_For_Ovito.xyz", "w");
	
	GenerateVerletLists(P, L, 0);

// ---------------------------------------------------------------------------------------- 

	// Main simulation loop
	for (int step = 0; step < sweeps; step++) {
		
		for (int i = 0; i < N; i++) {

			// Choses particle
			int index = (int)gsl_rng_uniform_int(generator, (unsigned long int)N);	

			// MC step
			if (gsl_rng_uniform(generator) < 0.5) {
				counter_acc_translation += MetropolisTranslation(index, P, L);
			} else {				
				counter_acc_rotation += MetropolisRotation(index, P, L);
			}

			// Volume step
			if (gsl_rng_uniform(generator) < 1. / (double)(N)) {
				double L_temp = L;
				L = VolumeMCStep(L, P);
				if(fabs(L_temp/L - 1.) > 1e-10) {
					counter_acc_volume++;
				}
			}
		}

		// Verlet-lists
		if(CheckVerlet(P, L)) {
			GenerateVerletLists(P, L, step);
			counter_verlet++;
		}

		// To monitor the acceptance rates of the individual steps
		if ((step + 1) % (sweeps / 10) == 0) {
			acc_rate[counter_store_current_rate][0] = 20.  * (double)counter_acc_translation / ((double)sweeps * (double)N);
			acc_rate[counter_store_current_rate][1] = 20.  * (double)counter_acc_rotation / ((double)sweeps * (double)N);
			acc_rate[counter_store_current_rate][2] = 10.  * (double)counter_acc_volume / ((double)sweeps);
			counter_acc_translation = 0;
			counter_acc_rotation = 0;
			counter_acc_volume = 0;
			counter_store_current_rate++;
		}

		// Prints position files
		if ((step + 1) % rate == 0) {
			PrintOvito(position_file_ovito, P, step, L);
		}
	}

// ---------------------------------------------------------------------------------------- 

	// Returns acceptance rates in console
	printf("end V: %f\t end L: %f\n", pow(L, (double)dim), L);
	printf("Verlet-lists were updated %ld times.\n", counter_verlet);	
	for (int i = 0; i < 10; i++) {
		printf("Acceptance rates:   tran. %1.4f   rot. %1.4f   vol. %1.4f\n", acc_rate[i][0], acc_rate[i][1], acc_rate[i][2]);
	}
	
	gsl_rng_free(generator);
	fclose(position_file_ovito);
	FreeVerLists(P);
	
	// Obtains the calendar time at the end of the simulation.
    time_t now;
    time(&now);

	// Calculating the duration of the simulation
    long simTime = now - start;
    printf("Duration of the simulation: %ld s\n", simTime);
	
    // End of simulation
    return 0;
}

// ---------------------------------------------------------------------------------------- 
