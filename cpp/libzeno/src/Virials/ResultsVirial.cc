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

#include "Virials/ResultsVirial.h"

using namespace zeno;

ResultsVirial::
    ResultsVirial(int numThreads, double refIntegral)
    :refAverage(NULL),
    refOverlapAverage(NULL),
    targetAverage(NULL),
    targetOverlapAverage(NULL),
    refAverageReduced(0),
    refOverlapAverageReduced(0),
    targetAverageReduced(0),
    targetOverlapAverageReduced(0),
    numThreads(numThreads),
    reduced(true),
    refNumSteps(NULL),
    targetNumSteps(NULL),
    refNumStepsReduced(0),
    targetNumStepsReduced(0),
    refIntegral(refIntegral),
    virialCoefficient(NULL),
    virialCoefficientReduced(0){
    refAverage = new Uncertain<double>[numThreads];
    refOverlapAverage = new Uncertain<double>[numThreads];
    targetAverage = new Uncertain<double>[numThreads];
    targetOverlapAverage = new Uncertain<double>[numThreads];
    refNumSteps = new long long[numThreads];
    targetNumSteps = new long long[numThreads];
    virialCoefficient = new Uncertain<double>[numThreads];

    for (int threadNum = 0; threadNum < numThreads; threadNum++) {
        refAverage[threadNum] = 0.0;
        refOverlapAverage[threadNum] = 0.0;
        targetAverage[threadNum] = 0.0;
        targetOverlapAverage[threadNum] = 0.0;
        refNumSteps[threadNum] = 0;
        targetNumSteps[threadNum] = 0;
        virialCoefficient[threadNum] = 0;
    }
}

ResultsVirial::
~ResultsVirial() {
    delete[] refAverage;
    delete[] refOverlapAverage;
    delete[] targetAverage;
    delete[] targetOverlapAverage;
    delete[] refNumSteps;
    delete[] targetNumSteps;
    delete[] virialCoefficient;
}

void
ResultsVirial::
putData(int threadNum,
        MeterOverlap<double> * refMeter,
        MeterOverlap<double> * targetMeter){

    reduced = false;

    double ** stats = refMeter->getStatistics();
    double ** blockCor = refMeter->getBlockCorrelation();
    refAverage[threadNum]=Uncertain<double>(stats[0][MeterOverlap<double>::AVG_AVG], squared(stats[0][MeterOverlap<double>::AVG_ERR]));
    refOverlapAverage[threadNum]=Uncertain<double>(stats[1][MeterOverlap<double>::AVG_AVG], squared(stats[1][MeterOverlap<double>::AVG_ERR]));
    double blockCov = blockCor[0][1] * refAverage[threadNum].getStdDev() * refOverlapAverage[threadNum].getStdDev();
    Uncertain<double>::setCovariance(refAverage[threadNum], refOverlapAverage[threadNum], blockCov);
    refNumSteps[threadNum] = refMeter->getNumSamples();
    stats = targetMeter->getStatistics();
    targetAverage[threadNum]=Uncertain<double>(stats[0][MeterOverlap<double>::AVG_AVG], squared(stats[0][MeterOverlap<double>::AVG_ERR]));
    targetOverlapAverage[threadNum]=Uncertain<double>(stats[1][MeterOverlap<double>::AVG_AVG], squared(stats[1][MeterOverlap<double>::AVG_ERR]))/refMeter->getAlpha()[0];
    blockCor = targetMeter->getBlockCorrelation();
    blockCov = blockCor[0][1] * targetAverage[threadNum].getStdDev() * targetOverlapAverage[threadNum].getStdDev();
    Uncertain<double>::setCovariance(targetAverage[threadNum], targetOverlapAverage[threadNum], blockCov);
    targetNumSteps[threadNum] = targetMeter->getNumSamples();
}

void
ResultsVirial::
reduce() {
    if (reduced) {
        return;
    }
    refAverageReduced = 0.0;
    refOverlapAverageReduced = 0.0;
    targetAverageReduced = 0.0;
    targetOverlapAverageReduced = 0.0;
    refNumStepsReduced = 0;
    targetNumStepsReduced = 0;
    virialCoefficientReduced = 0;

    for (int threadNum = 0; threadNum < numThreads; threadNum++) {
        refAverageReduced += refAverage[threadNum];
        refOverlapAverageReduced += refOverlapAverage[threadNum];
        refNumStepsReduced += refNumSteps[threadNum];
        targetAverageReduced += targetAverage[threadNum];
        targetOverlapAverageReduced += targetOverlapAverage[threadNum];
        targetNumStepsReduced += targetNumSteps[threadNum];
        virialCoefficientReduced += virialCoefficient[threadNum]/(double)numThreads;
    }
    reduced = true;
}

long long
ResultsVirial::
getNumSteps() const{
    assert(reduced);
    return refNumStepsReduced + targetNumStepsReduced;
}

double
ResultsVirial::
getRefFrac() const{
    assert(reduced);
    return (double)refNumStepsReduced/(refNumStepsReduced + targetNumStepsReduced);
}

Uncertain<double>
ResultsVirial::
getRefAverageReduced() const {
    assert(reduced);
    return refAverageReduced;
}

Uncertain<double>
ResultsVirial::
getRefOverlapAverageReduced() const {
    assert(reduced);
    return refOverlapAverageReduced;
}

Uncertain<double>
ResultsVirial::
getTargetAverageReduced() const {
    assert(reduced);
    return targetAverageReduced;
}

Uncertain<double>
ResultsVirial::
getTargetOverlapAverageReduced() const {
    assert(reduced);
    return targetOverlapAverageReduced;
}

double
ResultsVirial::
getRefIntegral() const {
    return refIntegral;
}

void
ResultsVirial::
putVirialCoefficient(int threadNum, double coefficient, double uncertainty){
    virialCoefficient[threadNum] = Uncertain<double>(coefficient, uncertainty);
}

Uncertain<double>
ResultsVirial::
getVirialCoefficientReduced() const {
    assert(reduced);
    return virialCoefficientReduced;
}
