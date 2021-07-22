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

#ifndef MC_MOVE_H
#define MC_MOVE_H

#include "ClusterSum.h"

namespace zeno {

/// Performs a monte carlo trial.
///
template <class T,
        class RandomNumberGenerator>
class IntegratorMSMC;

template <class T,
        class RandomNumberGenerator>
class MCMove {
 public:
    MCMove(IntegratorMSMC<T, RandomNumberGenerator> & integratorMSMC, ClusterSum<T> * clusterSum);
    virtual ~MCMove();
    virtual double doTrial(double oldValue, bool & accepted) = 0;
    double getStepSize();
    void setStepSize(double sS);

protected:
    IntegratorMSMC<T, RandomNumberGenerator> & integratorMSMC;
    ClusterSum<T> * clusterSum;
    double stepSize;
    void adjustStepSize();
    double maxStepSize;
    long numTrials, numAccepted;
    double chiSum;
    double adjustInterval;
    int lastAdjust;
    double adjustStep, minAdjustStep;
    bool verboseAdjust, tunable;
    
};

/// Sub class of MCMove to perform a monte carlo trial for translation.
///
template <class T,
        class RandomNumberGenerator>
class MCMoveTranslate : public MCMove<T, RandomNumberGenerator> {
public:
    MCMoveTranslate(IntegratorMSMC<T, RandomNumberGenerator> & integratorMSMC, ClusterSum<T> * clusterSum);
    ~MCMoveTranslate();

    double doTrial(double oldValue, bool & accepted);
 };

/// Sub class of MCMove to perform a monte carlo trial for rotation.
///
template <class T,
        class RandomNumberGenerator>
class MCMoveRotate : public MCMove<T, RandomNumberGenerator> {
public:
    MCMoveRotate(IntegratorMSMC<T, RandomNumberGenerator> & integratorMSMC, ClusterSum<T> * clusterSum);
    ~MCMoveRotate();

    double doTrial(double oldValue, bool & accepted);
};

/// Sub class of MCMove to perform a monte carlo trial for chain move.
///
template <class T,
        class RandomNumberGenerator>
class MCMoveChainVirial : public MCMove<T, RandomNumberGenerator> {
public:
    MCMoveChainVirial(IntegratorMSMC<T, RandomNumberGenerator> & integratorMSMC, ClusterSum<T> * clusterSum, double sigma);
    ~MCMoveChainVirial();

    double doTrial(double oldValue, bool & accepted);
protected:
    double sigma;
};

}
#endif

