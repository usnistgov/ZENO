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

#ifndef RESULTS_VIRIAL_H
#define RESULTS_VIRIAL_H

#include "Uncertain.h"
#include "Virials/MeterOverlap.h"

// ================================================================

namespace zeno {

/// Collects results from the virial coefficient computation.
///

class ResultsVirial{
public:
    ResultsVirial(int numThreads, double refIntegral);

    ~ResultsVirial();

    void putData(int threadNum,
                 MeterOverlap<double> * refMeter,
                 MeterOverlap<double> * targetMeter);

    void reduce();

    long long getNumSteps() const;
    double getRefIntegral() {
      return refIntegral;
    }
    double getRefFrac() const;
    Uncertain<double> getRefAverageReduced() const;
    Uncertain<double> getRefOverlapAverageReduced() const;
    Uncertain<double> getTargetAverageReduced() const;
    Uncertain<double> getTargetOverlapAverageReduced() const;
    double getRefIntegral() const;
    void putVirialCoefficient(int threadNum,
                              double coefficient,
                              double uncertainty);
    Uncertain<double> getVirialCoefficientReduced() const;
    static double squared(double x) {
      return x*x;
    }


private:
    Uncertain<double> * refAverage;
    Uncertain<double> * refOverlapAverage;
    Uncertain<double> * targetAverage;
    Uncertain<double> * targetOverlapAverage;
    Uncertain<double> refAverageReduced;
    Uncertain<double> refOverlapAverageReduced;
    Uncertain<double> targetAverageReduced;
    Uncertain<double> targetOverlapAverageReduced;
    int numThreads;
    bool reduced;
    long long * refNumSteps;
    long long * targetNumSteps;
    long long refNumStepsReduced;
    long long targetNumStepsReduced;
    const double refIntegral;
    Uncertain<double> * virialCoefficient;
    Uncertain<double> virialCoefficientReduced;
};

}
#endif //RESULTS_VIRIAL_H
