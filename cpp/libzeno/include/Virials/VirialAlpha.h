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

#ifndef VIRIAL_ALPHA_H
#define VIRIAL_ALPHA_H

#include "IntegratorMSMC.h"

// perhaps just take Clusters and make Meter and Average internally

template <class T,
        class RandomNumberGenerator>
class IntegratorMSMC;

template <class T,
        class RandomNumberGenerator>
class VirialAlpha {
protected:
    long stepCount, nextCheck;
    IntegratorMSMC<T, RandomNumberGenerator> & refIntegrator, & targetIntegrator;
    MeterOverlap<T> refMeter, targetMeter;
    double newAlpha, newAlphaErr, alphaCor, alphaSpan;
    bool allDone, verbose;
    double alphaStats[4];
public:
    VirialAlpha(IntegratorMSMC<T, RandomNumberGenerator> & refIntegrator,
                IntegratorMSMC<T, RandomNumberGenerator> & targetIntegrator,
                ClusterSum<T> & refClusterRef,
                ClusterSum<T> & refClusterTarget,
                ClusterSum<T> & targetClusterRef,
                ClusterSum<T> & targetClusterTarget);
    ~VirialAlpha();
    void setVerbose(bool newVerbose);
    void getNewAlpha(double & newAlpha, double & newAlphaErr, double & alphaCor);
    double * getAlphaStatistics();
    void setAlpha(double alphaCenter, double alphaSpan);
    void analyze(double & jBest);
    void runSteps(int steps);
    void run();
    bool getAllDone();
};


#endif //VIRIAL_ALPHA_H
