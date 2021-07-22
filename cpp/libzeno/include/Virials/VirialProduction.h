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
#ifndef VIRIAL_PRODUCTION_H
#define VIRIAL_PRODUCTION_H

#include "IntegratorMSMC.h"

// perhaps just take Clusters and make Meter and Average internally

template <class T,
        class RandomNumberGenerator>
class IntegratorMSMC;

template <class T,
        class RandomNumberGenerator>
class VirialProduction {
protected:
    IntegratorMSMC<T, RandomNumberGenerator> & refIntegrator, & targetIntegrator;
    MeterOverlap<T> refMeter, targetMeter;
    double idealTargetFraction;
    double **refStats, **refBCStats, **refRatioStats;
    double **targetStats, **targetBCStats, **targetRatioStats;
    double **fullBCStats;
    double alphaStats[2];
    double **fullStats;
    double refIntegral;
    bool disposed;
    long refSteps, targetSteps;
public:
    VirialProduction(IntegratorMSMC<T, RandomNumberGenerator> &refIntegrator,
                     IntegratorMSMC<T, RandomNumberGenerator> &targetIntegrator,
                     ClusterSum<T> &refClusterRef,
                     ClusterSum<T> &refClusterTarget,
                     ClusterSum<T> &targetClusterRef,
                     ClusterSum<T> &targetClusterTarget,
                     double alpha, double refIntegral);
    ~VirialProduction();
    void analyze();
    void printResults(const char **targetNames);
    void getResults();
    void runSteps(long long numSteps);
    double** getFullStats();
    double* getAlphaStats();
    double** getRefStats();
    double** getTargetStats();
    double** getRefBCStats();
    double** getTargetBCStats();
    double** getRefRatioStats();
    double** getTargetRatioStats();
    double** getFullBCStats();
    MeterOverlap<T> * getRefMeter();
    MeterOverlap<T> * getTargetMeter();
};

#endif //VIRIAL_PRODUCTION_H
