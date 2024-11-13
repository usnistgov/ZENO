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
    ResultsVirial(int numThreads, int numValues, double refIntegral);

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
    Uncertain<double> getTargetAverageReduced(int iValue) const;
    Uncertain<double> getOverlapRatioAverageReduced() const;
    double getRefIntegral() const;
    void putOverlapRatio(int threadNum,
                         double overlapRatio,
                         double uncertainty);
    void putVirialCoefficient(int threadNum,
                              int iValue,
                              double coefficient,
                              double uncertainty);
    Uncertain<double> getVirialCoefficientReduced(int iValue) const;
    static double squared(double x) {
      return x*x;
    }
    int getNumValues() const {
      return virialCoefficientReduced.size();
    }

    void setOrder(int order) {
      this->order = order;
    }
    int getOrder() const {
      return this->order;
    }

private:
    Uncertain<double> * refAverage;
    Uncertain<double> * targetAverage;
    Uncertain<double> * overlapRatioAverage;
    Uncertain<double> refAverageReduced;
    Uncertain<double> targetAverageReduced;
    Uncertain<double> overlapRatioAverageReduced;
    int order;
    int numThreads;
    bool reduced;
    long long * refNumSteps;
    long long * targetNumSteps;
    long long refNumStepsReduced;
    long long targetNumStepsReduced;
    const double refIntegral;
    std::vector<Uncertain<double>> * virialCoefficient;
    std::vector<Uncertain<double>> virialCoefficientReduced;
};

}
#endif //RESULTS_VIRIAL_H
