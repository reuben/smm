#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdio.h>
#include <math.h>
#include "smm.h"

/* Compile with gcc -std=c99 test.c smm.c -lgsl -lm */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define SIZE 1000   // max number of sensors.

typedef struct { double min, avg, max; } Pair;

Pair stats(double v[], int n)
{
    Pair p;
    p.min = p.max = p.avg = v[0];
    for(int i=1; i<n; i++) {
        if (v[i] < p.min) p.min = v[i];
        if (v[i] > p.max) p.max = v[i];
        p.avg += v[i];
    }
    p.avg /= n;
    return p;
}

/* SmmTest:
   * Init the library for Meandering Jet
   * Deploy the sensors
   * Trace the fastest along the x axis...(the flow).
*/

int main(int argc, char *argv[])
{
    double x[SIZE];
    double y[SIZE];

    int p = atoi(argv[1]);

    int deployed = 0;

    gsl_rng *r = gsl_rng_alloc(gsl_rng_ranlxs1);
    int seed = 1214314;
    gsl_rng_set(r, seed);

    const double c = 0.12;
    const double k = 2.0*M_PI/7.5;
    const double T = (2*M_PI)/(c*k);
    double randomshift = gsl_rng_uniform(r)*T;

    // printf("# seed = %d, phase = %f\n", seed, randomshift);

    const int Streamfunction = SMM_MNDRJET;
    double params[SMM_MXPRMS]  = { 1.2, 0.12, 2.0*M_PI/7.5, 0.4, 0.3, randomshift };

    // Passo di Integrazione
    const double SmmTimeStep = 0.01;
    // 1 t = 0.03 Days.
    const double DimensionalTime = 0.03;
    // Durata in secondi del singolo passo = 25.92
    const double MicroStepTime = SmmTimeStep*24*60*60 * DimensionalTime;

    SMM_init(SIZE, Streamfunction, params, SMM_RK2, SmmTimeStep);

    // initial deployment 200 sensors in the domain [0,4] x [-2,2] km

    int n = atoi(argv[2]);
    for (int i=0; i<n; i++) {
        x[i] = gsl_rng_uniform(r) * 4.0;
        y[i] = gsl_rng_uniform(r) * 4.0 - 2.0;
        printf("$node_(%d) set X_ %lf\n$node_(%d) set Y_ %lf\n", i, x[i]*1000, i, y[i]*1000);
    }
    deployed += n;
    SMM_deploy_nodes(n);

    // run for one day -> 24 h x 60x60 s.
    // trace the fastest particle;

    // printf("# Time - min(X) - max(X)\n");
    // Pair xp = stats(x,n);
    // printf("%7.2f\t %.2f\t %.2f  %.2f\n", 0.0,
    //         xp.min*1000, xp.avg*1000, xp.max*1000);

    double time;
    for (int i=0; i<ceil((p*60*60)/25.92); i++) {
        SMM_Move_Sensors(x,y);
        // Pair xp = stats(x,n);
        // time in seconds..
        time = SMM_get_current_time() * (MicroStepTime/SmmTimeStep);
        for (int i = 0; i < n; ++i) {
            printf("$ns_ at %lf \"$node_(%d) setdest %lf %lf 0.0\"\n", time, i, x[i]*1000, y[i]*1000);
        }
        // printf("%.2f,%.2f,%.2f\n", time, x[p], y[p]);
        // printf("%7.2f\t %.2f\t %.2f  %.2f\n",
        //         time, xp.min*1000, xp.avg*1000, xp.max*1000);
    }
}
