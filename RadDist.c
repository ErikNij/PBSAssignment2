#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "RadDist.h"
#include "constants.h"
#include "structs.h"
#include "nbrlist.h"

void radDist(struct Parameters *p_parameters,struct Nbrlist *p_nbrlist,FILE *fpt)
{
    struct Pair *nbr = p_nbrlist->nbr;
    const size_t num_nbrs = p_nbrlist->num_nbrs;
    int *binRes;
    double *binResCorr;

    int numbBins = 100;
    float binMin = 0.0;
    float binMax = 10.0;
    float dbin = (binMax - binMin)/numbBins;
    struct DeltaR rij;
    int binNumb;
    double shellVol;
    double simVol = p_parameters->L.x * p_parameters->L.y * p_parameters->L.z;
    double expectedParts = p_parameters->num_part / simVol;

    binRes = (int *)malloc(numbBins * sizeof(int));
    binResCorr = (double *)malloc(numbBins * sizeof(double));
    
    for (int i = 0; i < numbBins; i++)
    {
        binRes[i] = 0;
        binResCorr[i] = 0;
    }

    for (size_t k = 0; k < num_nbrs; k++)
    {
        rij = nbr[k].rij;
        binNumb = floor(sqrt(rij.sq) / dbin);
        if (binNumb > numbBins)
        {
            binNumb = numbBins;
        }
        binRes[binNumb]++;
    }

    for (int i = 0; i < numbBins; i++)
    {
        shellVol =  (4/3*PI*pow((i+1)*dbin,3)-(4/3*PI*pow(i*dbin,3)));
        binResCorr[i] = binRes[i]/shellVol/expectedParts;
        fprintf(fpt,"%lf,",binResCorr[i]);
    }
    fprintf(fpt,"\n");



}