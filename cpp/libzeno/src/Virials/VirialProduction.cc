// ================================================================
//
/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
//
// ================================================================

// ================================================================
//
// Authors:
// Created:
//
// ================================================================

#include "Virials/VirialProduction.h"
#include "Virials/alloc2D.h"

// perhaps just take Clusters and make Meter and Average internally

template <class T,
        class RandomNumberGenerator>
VirialProduction<T, RandomNumberGenerator>::
VirialProduction(IntegratorMSMC<T, RandomNumberGenerator> & rIntegrator,
                 IntegratorMSMC<T, RandomNumberGenerator> & tIntegrator,
                 ClusterSum<T> &refClusterRef,
                 ClusterSum<T> &refClusterTarget,
                 ClusterSum<T> &targetClusterRef,
                 ClusterSum<T> &targetClusterTarget,
                 double alpha, double ri) :
                 refIntegrator(rIntegrator), targetIntegrator(tIntegrator),
                 refMeter(MeterOverlap<T>(refClusterTarget, alpha, 0, 1)),
                 targetMeter(MeterOverlap<T>(targetClusterRef, 1/alpha, 0, 1)),
                 idealTargetFraction(0.5), refIntegral(ri), refSteps(0), targetSteps(0) {
    int numTargets = targetMeter.getNumData();
    fullStats = (double**)malloc2D(numTargets-1, 2, sizeof(double));
    fullBCStats = (double**)malloc2D(numTargets-1, numTargets-1, sizeof(double));
    refIntegrator.setMeter(&refMeter);
    targetIntegrator.setMeter(&targetMeter);
}

template <class T,
        class RandomNumberGenerator>
VirialProduction<T, RandomNumberGenerator>::
  ~VirialProduction() {
  free2D((void**)fullStats);
  free2D((void**)fullBCStats);
}

template <class T,
        class RandomNumberGenerator>
void VirialProduction<T, RandomNumberGenerator>::
analyze() {
    refStats = refMeter.getStatistics();
    targetStats = targetMeter.getStatistics();
    if (std::isnan(targetStats[0][MeterOverlap<T>::AVG_ERR])) {
        idealTargetFraction = 1;
        return;
    }
    if (std::isnan(refStats[0][MeterOverlap<T>::AVG_ERR])) {
        idealTargetFraction = 0;
        return;
    }
    int numTargets = targetMeter.getNumData();
    double alpha = refMeter.getAlpha()[0];
    alphaStats[0] = refStats[1][MeterOverlap<T>::AVG_AVG]/targetStats[numTargets-1][MeterOverlap<T>::AVG_AVG]*alpha;
    alphaStats[1] = alpha*MeterOverlap<T>::ratioErr(refStats[1][MeterOverlap<T>::AVG_AVG], refStats[1][MeterOverlap<T>::AVG_ERR],
                                                 targetStats[numTargets-1][MeterOverlap<T>::AVG_AVG], targetStats[numTargets-1][MeterOverlap<T>::AVG_ERR], 0);

    refRatioStats = refMeter.getRatioStatistics();
    targetRatioStats = targetMeter.getRatioStatistics();
    refBCStats = refMeter.getBlockCorrelation();
    double **cij = targetMeter.getRatioCorrelation();
    targetBCStats = targetMeter.getBlockCorrelation();
    double vd = refRatioStats[0][MeterOverlap<T>::AVG_AVG];
    double ed = refRatioStats[0][MeterOverlap<T>::AVG_ERR];
    for (int i=0; i<numTargets-1; i++) {
        fullStats[i][0] = refIntegral * targetRatioStats[i][MeterOverlap<T>::AVG_AVG] / vd * alpha;
        fullStats[i][1] = fabs(refIntegral) * MeterOverlap<T>::ratioErr(targetRatioStats[i][MeterOverlap<T>::AVG_AVG], targetRatioStats[i][MeterOverlap<T>::AVG_ERR], vd, ed,0) * alpha;
        double vi = targetRatioStats[i][MeterOverlap<T>::AVG_AVG];
        double ei = targetRatioStats[i][MeterOverlap<T>::AVG_ERR];
        for (int j=0; j<numTargets-1; j++) {
            if (j==i) {
                fullBCStats[i][i] = 1;
                continue;
            }
            double vj = targetRatioStats[j][MeterOverlap<T>::AVG_AVG];
            double ej = targetRatioStats[j][MeterOverlap<T>::AVG_ERR];
            fullBCStats[i][j] = fullBCStats[j][i] = MeterOverlap<T>::ratioCor(vi, vj, vd, ei, ej, ed, cij[i][j], 0, 0);
        }
    }

    double oldFrac = ((double)targetSteps)/(targetSteps + refSteps);
    double refErrorRatio = refRatioStats[0][MeterOverlap<T>::AVG_ERR]/fabs(refRatioStats[0][MeterOverlap<T>::AVG_AVG]);
    if (std::isnan(refErrorRatio) || refErrorRatio > 1) refErrorRatio = 1;
    double targetErrorRatio = targetRatioStats[0][MeterOverlap<T>::AVG_ERR]/fabs(targetRatioStats[0][MeterOverlap<T>::AVG_AVG]);
    if (std::isnan(targetErrorRatio) || targetErrorRatio > 1) targetErrorRatio = 1;
    idealTargetFraction = 1.0/(1 + refErrorRatio/targetErrorRatio * sqrt((1-oldFrac)/oldFrac));
}

template <class T,
        class RandomNumberGenerator>
void VirialProduction<T, RandomNumberGenerator>::
printResults(const char **targetNames) {
    int numTargets = targetMeter.getNumData();
    printf("final reference step fraction: %5.4f\n", 1-idealTargetFraction);
    printf("actual reference step fraction: %5.4f\n", ((double)refSteps)/(refSteps+targetSteps));
    printf("reference blocks: %lld of size %lld\n", refMeter.getBlockCount(), refMeter.getBlockSize());
    printf("target blocks: %lld of size %lld\n", targetMeter.getBlockCount(), targetMeter.getBlockSize());
    printf("alpha check:               % 22.15e  error: %12.5e\n", alphaStats[0], alphaStats[1]);
    printf("full average:              % 22.15e  error: %12.5e\n", fullStats[0][0], fullStats[0][1]);
    printf("reference ratio:           % 22.15e  error: %12.5e   cor: % 7.5f\n", refRatioStats[0][MeterOverlap<T>::AVG_AVG], refRatioStats[0][MeterOverlap<T>::AVG_ERR], refBCStats[0][1]);
    printf("reference average:         % 22.15e  error: %12.5e  acor: % 7.5f\n", refStats[0][MeterOverlap<T>::AVG_AVG], refStats[0][MeterOverlap<T>::AVG_ERR], refStats[0][MeterOverlap<T>::AVG_ACOR]);
    printf("reference overlap average: % 22.15e  error: %12.5e  acor: % 7.5f\n", refStats[1][MeterOverlap<T>::AVG_AVG], refStats[1][MeterOverlap<T>::AVG_ERR], refStats[1][MeterOverlap<T>::AVG_ACOR]);
    printf("target ratio:              % 22.15e  error: %12.5e   cor: % 7.5f\n", targetRatioStats[0][MeterOverlap<T>::AVG_AVG], targetRatioStats[0][MeterOverlap<T>::AVG_ERR], targetBCStats[0][numTargets-1]);
    printf("target average:            % 22.15e  error: %12.5e  acor: % 7.5f\n", targetStats[0][MeterOverlap<T>::AVG_AVG], targetStats[0][MeterOverlap<T>::AVG_ERR], targetStats[0][MeterOverlap<T>::AVG_ACOR]);
    printf("target overlap average:    % 22.15e  error: %12.5e  acor: % 7.5f\n", targetStats[numTargets-1][MeterOverlap<T>::AVG_AVG], targetStats[numTargets-1][MeterOverlap<T>::AVG_ERR], targetStats[numTargets-1][MeterOverlap<T>::AVG_ACOR]);
    if (numTargets == 2) return;
    for (int i=1; i<numTargets-1; i++) {
        char name[40];
        if (targetNames && strlen(targetNames[i]) > 11) {
            fprintf(stderr, "truncating name %s to 11 characters\n", targetNames[i]);
        }
        if (targetNames && targetNames[i]) snprintf(name, 39, "%.11s average:", targetNames[i]);
        else snprintf(name, 39, "extra %d average:", i);
        printf("%-26s % 22.15e  error: %12.5e  acor: % 7.5f\n", name, targetStats[i][MeterOverlap<T>::AVG_AVG], targetStats[i][MeterOverlap<T>::AVG_ERR], targetStats[i][MeterOverlap<T>::AVG_ACOR]);
        if (targetNames && targetNames[i]) snprintf(name, 39, "%.11s ratio average:", targetNames[i]);
        else snprintf(name, 39, "extra %d ratio average:", i);
        printf("%-26s % 22.15e  error: %12.5e   cor: % 7.5f\n", name, targetRatioStats[i][MeterOverlap<T>::AVG_AVG], targetRatioStats[i][MeterOverlap<T>::AVG_ERR], targetBCStats[i][numTargets-1]);
        if (targetNames && targetNames[i]) snprintf(name, 39, "full %.11s average:", targetNames[i]);
        else snprintf(name, 39, "full extra %d average:", i);
        printf("%-26s % 22.15e  error: %12.5e\n", name, fullStats[i][0], fullStats[i][1]);
    }
    if (numTargets == 1) return;
    printf("Target Correlation:\n");
    for (int i=0; i<numTargets-1; i++) {
        for (int j=0; j<numTargets-1; j++) {
            printf(" % 8.5f", targetBCStats[i][j]);
        }
        printf("\n");
    }
    printf("Full Correlation:\n");
    for (int i=0; i<numTargets-1; i++) {
        for (int j=0; j<numTargets-1; j++) {
            printf(" % 8.5f", fullBCStats[i][j]);
        }
        printf("\n");
    }
}

template <class T,
        class RandomNumberGenerator>
double** VirialProduction<T, RandomNumberGenerator>::
getFullStats() {
    return fullStats;
}

template <class T,
        class RandomNumberGenerator>
double* VirialProduction<T, RandomNumberGenerator>::
getAlphaStats() {
    return alphaStats;
}

template <class T,
        class RandomNumberGenerator>
double** VirialProduction<T, RandomNumberGenerator>::
getRefStats() {
    return refStats;
}

template <class T,
        class RandomNumberGenerator>
double** VirialProduction<T, RandomNumberGenerator>::
getTargetStats() {
    return targetStats;
}

template <class T,
        class RandomNumberGenerator>
double** VirialProduction<T, RandomNumberGenerator>::
getRefBCStats() {
    return refBCStats;
}

template <class T,
        class RandomNumberGenerator>
double** VirialProduction<T, RandomNumberGenerator>::
getTargetBCStats() {
    return targetBCStats;
}

template <class T,
        class RandomNumberGenerator>
double** VirialProduction<T, RandomNumberGenerator>::
getRefRatioStats() {
    return refRatioStats;
}

template <class T,
        class RandomNumberGenerator>
double** VirialProduction<T, RandomNumberGenerator>::
getTargetRatioStats() {
    return targetRatioStats;
}

template <class T,
        class RandomNumberGenerator>
double** VirialProduction<T, RandomNumberGenerator>::
getFullBCStats() {
    return fullBCStats;
}

template <class T,
        class RandomNumberGenerator>
void VirialProduction<T, RandomNumberGenerator>::
runSteps(long long numSteps) {
    long totalSteps = refSteps + targetSteps;
    long subSteps = 100 + totalSteps/1000;
    if (subSteps > numSteps) subSteps = numSteps;
    long thisSteps = 0;
    while (thisSteps < numSteps) {
        bool runRef = true;
        if (totalSteps > 0) {
            double tFrac = ((double)targetSteps)/totalSteps;
            double idf = std::max(std::min(idealTargetFraction,0.99),0.01);
            runRef = tFrac > idf;
        }

        if (runRef) {
            refIntegrator.doStep(subSteps);
            refSteps += subSteps;
        }
        else {
            targetIntegrator.doStep(subSteps);
            targetSteps += subSteps;
        }

        analyze();
        totalSteps += subSteps;
        thisSteps += subSteps;

        subSteps = 100 + totalSteps/1000;
        if (subSteps > numSteps-thisSteps) subSteps = numSteps - thisSteps;
    }
}

template <class T,
        class RandomNumberGenerator>
MeterOverlap<T> *
VirialProduction<T, RandomNumberGenerator>::
getRefMeter() {
    return &refMeter;
}

template <class T,
        class RandomNumberGenerator>
MeterOverlap<T> *
VirialProduction<T, RandomNumberGenerator>::
getTargetMeter() {
    return &targetMeter;
}
