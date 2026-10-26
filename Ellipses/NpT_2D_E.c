#include "NpT_2D_E.h"
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
int dim = DIM;												// Dimension of the system. Should be fixed to 2
int N = NP;													// Number of particles
double p = PRESSURE;										// Pressure (set by user)
double T = TEMP;											// Temperature (should be set to 1, not tested otherwise)
int sweeps = SWEEPS;										// Speeps or MC steps
int rate = RATE;											// Rate which dictates when data is written
double a = AA;												// Longer semi-axis
double b = BB;												// Shorter semi-axis
long seed = SEED;											// Random seed
double delta = DELTA;										// Displacement translation
double delta_r = DELTAR;									// Displacement rotation
double delta_V = DELTAV;									// Displacement volume
double r_v = RV;											// Radius for Verlet lists

// ----------------------------------------------------------------------------------------
// Main
// ----------------------------------------------------------------------------------------

int main(int argc, char *argv[]) {
	
	// For the calculation of the total simulation time. 
    time_t start;
    time(&start);
	
	PrintWelcome();
	
	// The program only runs with "inputfile" in the same folder. 
    // A new inputfile can be generated with "./NpT_2D_E.exe Setup"
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
	ReadInput();	// See inputdata.c
	PrintInput();	// See printdisplay.c
	CheckForConflict();	// See inputdata.c

	// Sets initial system area. V_init_min is chosen through tests and can be varied if needed.
	double V = (double)N * T / p;
	double V_init_min = 4.7 * (double)N * T / 3.;
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
	InitRandom(P, L);	// See init.c
	
	double acc_rate[10][3];
	long counter_verlet = 1;
	long counter_acc_translation = 0;
	long counter_acc_rotation = 0;
	long counter_acc_volume = 0;
	long counter_store_current_rate = 0;
	
	FILE *position_file_ovito = NULL;
    position_file_ovito = fopen("Position_For_Ovito.xyz", "w");
	
	GenerateVerletLists(P, L, 0);	// See verletlist.c

	// ---------------------------------------------------------------------------------------- 

	// Main simulation loop
	for (int step = 0; step < sweeps; step++) {
		
		for (int i = 0; i < N; i++) {

			// Choses particle
			int index = (int)gsl_rng_uniform_int(generator, (unsigned long int)N);	

			// MC step; see mc_steps.c
			if (gsl_rng_uniform(generator) < 0.5) {
				counter_acc_translation += MetropolisTranslation(index, P, L);
			} else {				
				counter_acc_rotation += MetropolisRotation(index, P, L);
			}

			// Volume step
			if (gsl_rng_uniform(generator) < 1. / (double)(N)) {
				double L_temp = L;
				L = VolumeMCStep(L, P);	// See mc_steps.c
				if(fabs(L_temp/L - 1.) > 1e-10) {
					counter_acc_volume++;
				}
			}
		}

		// Verlet-lists; see verletlists.c
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

		// Prints position files; see printdata.c
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
	FreeVerLists(P); // See verletlists.c
	
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
