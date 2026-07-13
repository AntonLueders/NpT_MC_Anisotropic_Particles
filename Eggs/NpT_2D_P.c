#include "NpT_2D_P.h"
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
double b = BB;
long seed = SEED;
double delta = DELTA;
double delta_r = DELTAR;
double delta_V = DELTAV;
double r_v = RV;


int main(int argc, char *argv[]) {
	
	// For the calculation of the total simulation time. 
    time_t start;
    time(&start);
	
	PrintWelcome();
	
	// The program only runs with "inputfile" in the same folder. 
    // A new inputfile can be generated with "./NpT_SC.exe Setup"
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
	
	ReadInput();
	PrintInput();
	CheckForConflict();
	
	double V = (double)N * T / p;
	double V_init_min = 5.5 * (double)N * T / 3.;
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
	
	/*P[0].r[0] = 0;
	P[0].r[1] = 0;
	
	P[1].r[0] = -1.3;
	P[1].r[1] = 0.2;
	
	P[0].e[0] =  -0.5596;    
	P[0].e[1] =  -0.8288;
	
	P[1].e[0] =  -0.4752;
	P[1].e[1] =  -0.8799;
	
	int test[2];
	IdentifyEnds(&P[0], &P[1], L, test);
	printf("%d\t%d\n", test[0], test[1]);*/
	
	for (int step = 0; step < sweeps; step++) {
		
		for (int i = 0; i < N; i++) {
			
			int index = (int)gsl_rng_uniform_int(generator, (unsigned long int)N);	
			
			if (gsl_rng_uniform(generator) < 0.5) {
				counter_acc_translation += MetropolisTranslation(index, P, L);
			} else {				
				counter_acc_rotation += MetropolisRotation(index, P, L);
			}

			if (gsl_rng_uniform(generator) < 1. / (double)(N)) {
				double L_temp = L;
				
				L = VolumeMCStep(L, P);
				if(fabs(L_temp/L - 1.) > 1e-10) {
					counter_acc_volume++;
				}
				
			}
		}
		
		if(CheckVerlet(P, L)) {
			GenerateVerletLists(P, L, step);
			counter_verlet++;
		}
		
		if ((step + 1) % (sweeps / 10) == 0) {
			acc_rate[counter_store_current_rate][0] = 20.  * (double)counter_acc_translation / ((double)sweeps * (double)N);
			acc_rate[counter_store_current_rate][1] = 20.  * (double)counter_acc_rotation / ((double)sweeps * (double)N);
			acc_rate[counter_store_current_rate][2] = 10.  * (double)counter_acc_volume / ((double)sweeps);
			counter_acc_translation = 0;
			counter_acc_rotation = 0;
			counter_acc_volume = 0;
			counter_store_current_rate++;
		}
		
		if ((step + 1) % rate == 0) {
			PrintOvito(position_file_ovito, P, step, L);
		}
	}
	
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
