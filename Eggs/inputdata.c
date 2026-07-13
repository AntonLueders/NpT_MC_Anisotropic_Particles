#include "inputdata.h"

// ----------------------------------------------------------------------------------------

// Prints inputfile if "Setup" is used.
void PrintInputfile() {

    printf("Write intput file using the default values. Please wait...\n");

    printf("\n");

    FILE *file = fopen("inputfile", "w");
    fprintf(file, "N\t%d\n", N);
    fprintf(file, "p\t%1.6f\n", p);
	fprintf(file, "T\t%1.6f\n", T);
	fprintf(file, "sweeps\t%d\n", sweeps);
	fprintf(file, "rate\t%d\n", rate);
	fprintf(file, "a\t%1.6f\n", a);
	fprintf(file, "b\t%1.6f\n", b);
	fprintf(file, "seed\t%ld\n", seed);
	fprintf(file, "delta\t%1.6f\n", delta);
	fprintf(file, "delta_r\t%1.6f\n", delta_r);
	fprintf(file, "delta_V\t%1.6f\n", delta_V);
	fprintf(file, "r_v\t%1.6f\n", r_v);

    fclose(file);
}

// ----------------------------------------------------------------------------------------

// Reads inputfile. First, the parameters are stored in the "temp" variables are generated. 
// If reading the input was successfully, the global variables are set to the temp values.
void ReadInput() {

    // If there is no inputfile, an error occures.
    FILE *input = fopen("inputfile", "r");
    if (input == NULL) {
        printf("ERROR: No input found. Use Setup option.\n");

        printf("\n");

        exit(1);
    }
	
	// Temp variables
    char dummy[40] = " ";
    char temp[40] = " ";
	int N_temp = N;
	double p_temp = p;
	double T_temp = T;
	int sweeps_temp = sweeps;
	int rate_temp = rate;
	double a_temp = a;
	double b_temp = b;
	long seed_temp = (long)seed;
	double delta_temp = delta;
	double delta_r_temp = delta_r;
	double delta_V_temp = delta_V;
	double r_v_temp = r_v;

	
	// Reading the inputfile
    if (fscanf(input, "%s\t%d%[^\n]\n", temp, &N_temp, dummy) == 2) {
        N = N_temp;
    } else {
        printf("ERROR: Check input or use Setup option. (N)\n");
        exit(1);
    }
	if (fscanf(input, "%s\t%lf%[^\n]\n", temp, &p_temp, dummy) == 2) {
        p = p_temp;
    } else {
        printf("ERROR: Check input or use Setup option. (p)\n");
        exit(1);
    }	
	if (fscanf(input, "%s\t%lf%[^\n]\n", temp, &T_temp, dummy) == 2) {
        T = T_temp;
    } else {
        printf("ERROR: Check input or use Setup option. (T)\n");
        exit(1);
    }
	if (fscanf(input, "%s\t%d%[^\n]\n", temp, &sweeps_temp, dummy) == 2) {
        sweeps = sweeps_temp;
    } else {
        printf("ERROR: Check input or use Setup option. (sweeps)\n");
        exit(1);
    }
	if (fscanf(input, "%s\t%d%[^\n]\n", temp, &rate_temp, dummy) == 2) {
        rate = rate_temp;
    } else {
        printf("ERROR: Check input or use Setup option. (rate)\n");
        exit(1);
    }
	if (fscanf(input, "%s\t%lf%[^\n]\n", temp, &a_temp, dummy) == 2) {
        a = a_temp;
    } else {
        printf("ERROR: Check input or use Setup option. (a)\n");
        exit(1);
    }
		if (fscanf(input, "%s\t%lf%[^\n]\n", temp, &b_temp, dummy) == 2) {
        b = b_temp;
    } else {
        printf("ERROR: Check input or use Setup option. (b)\n");
        exit(1);
    }
    if (fscanf(input, "%s\t%ld%[^\n]\n", temp, &seed_temp, dummy) == 2) {
        seed = (long)seed_temp;
    } else {
        printf("ERROR: Check input or use Setup option. (seed)\n");
        exit(1);
    }	
	if (fscanf(input, "%s\t%lf%[^\n]\n", temp, &delta_temp, dummy) == 2) {
        delta = delta_temp;
    } else {
        printf("ERROR: Check input or use Setup option. (delta)\n");
        exit(1);
    }
	if (fscanf(input, "%s\t%lf%[^\n]\n", temp, &delta_r_temp, dummy) == 2) {
        delta_r = delta_r_temp;
    } else {
        printf("ERROR: Check input or use Setup option. (delta_r)\n");
        exit(1);
    }
	if (fscanf(input, "%s\t%lf%[^\n]\n", temp, &delta_V_temp, dummy) == 2) {
        delta_V = delta_V_temp;
    } else {
        printf("ERROR: Check input or use Setup option. (delta_V)\n");
        exit(1);
    }
	if (fscanf(input, "%s\t%lf%[^\n]\n", temp, &r_v_temp, dummy) == 2) {
        r_v = r_v_temp;
    } else {
        printf("ERROR: Check input or use Setup option. (r_v)\n");
        exit(1);
    }	
	
	fclose(input);
}

// Checks for conflicting input
void CheckForConflict() {
    
	if (N < 1) {
        printf("ERROR: Use at least 1 particle.\n");
        exit(1);
    }
    if (p < 0.) {
        printf("ERROR: Negative pressure is not defined.\n");
        exit(1);
    }
    
    if (T < 0.) {
        printf("ERROR: Negative temperature is not defined.\n");
        exit(1);        
    }
    
    if (rate > sweeps) {
        printf("ERROR: Rate must be larger than sweeps.\n");
        exit(1);    
    }
    
    if (a <= b) {
        printf("ERROR: a must be larger than b.\n");
        exit(1);
    }
    
    if (r_v < 2. * a + 3. * delta) {
        printf("WARNING: Particles can move a larger distance than r_v.\n");
    }

    if (a / b > 20. && a / b < 30.) {
		printf("WARNING: The computations might be unstable for a /b > 20.\n");
		printf("WARNING: Please use a /b <= 20.\n");
	} else if (a / b >= 30.) {
		printf("ERROR: The computations are designed for a / b <= 20.\n");
        exit(1);  
	}
}
