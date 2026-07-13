#include "printdisplay.h"

// ----------------------------------------------------------------------------------------

// Welcome screen with the information of the program. (Is printed to the console).
void PrintWelcome() {

    printf("\n");

    printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
           "++++"
           "+++++++\n");
    printf("+                                                                  "
           "    "
           "      +\n");
    printf("+             NpT Monte Carlo simulations of hard convex bodies    "
           "    "
           "      +\n");
    printf("+             2D Spherocylinder version                            "
           "    "
           "      +\n");
    printf("+                                                                  "
           "    "
           "      +\n");
    printf("+             Author:   Anton Lueders                              "
           "    "
           "      +\n");
    printf("+             Date:     20.02.2025                                 "
           "    "
           "      +\n");
    printf("+                                                                  "
           "    "
           "      +\n");
    printf("+             AG Franosch                                          "
           "    "
           "      +\n");
    printf("+             BioNano-Physics                                      "
           "    "
           "      +\n");
    printf("+             Univesity of Innsbruck                               "
           "    "
           "      +\n");
    printf("+                                                                  "
           "    "
           "      +\n");
    printf("+             Last edited: 12.09.2025                              "
           "    "
           "      +\n");
    printf("+                                                                  "
           "    "
           "      +\n");
    printf("+             V2.0.4                                               "
           "          +\n");
    printf("+                                                                  "
           "    "
           "      +\n");
    printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
           "++++"
           "+++++++\n");

    printf("\n");
}

// ----------------------------------------------------------------------------------------

// Prints the input data (global variables) to the console.
void PrintInput() {

    printf("The simulation was started using the following input:\n");

    printf("\n");
	
	printf("N\t%d\n", N);
    printf("p\t%1.6f\n", p);
	printf("T\t%1.6f\n", T);
	printf("sweeps\t%d\n", sweeps);
	printf("rate\t%d\n", rate);
	printf("q\t%1.6f\n", q);
	printf("sigma\t%1.6f\n", sigma);
	printf("seed\t%ld\n", seed);
	printf("delta\t%1.6f\n", delta);
	printf("delta_r\t%1.6f\n", delta_r);
	printf("delta_V\t%1.6f\n", delta_V);
	printf("r_v\t%1.6f\n", r_v);
}

// ----------------------------------------------------------------------------------------
