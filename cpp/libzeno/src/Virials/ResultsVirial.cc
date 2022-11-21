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
  : refAverage(NULL),
    targetAverage(NULL),
    overlapRatioAverage(NULL),
    refAverageReduced(0),
    targetAverageReduced(0),
    overlapRatioAverageReduced(0),
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
    targetAverage = new Uncertain<double>[numThreads];
    overlapRatioAverage = new Uncertain<double>[numThreads];
    refNumSteps = new long long[numThreads];
    targetNumSteps = new long long[numThreads];
    virialCoefficient = new Uncertain<double>[numThreads];

    for (int threadNum = 0; threadNum < numThreads; threadNum++) {
        refAverage[threadNum] = 0.0;
        targetAverage[threadNum] = 0.0;
        overlapRatioAverage[threadNum] = 0.0;
        refNumSteps[threadNum] = 0;
        targetNumSteps[threadNum] = 0;
        virialCoefficient[threadNum] = 0;
    }
}

ResultsVirial::
~ResultsVirial() {
    delete[] refAverage;
    delete[] targetAverage;
    delete[] overlapRatioAverage;
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
    Uncertain<double> refOverlapAverage = Uncertain<double>(stats[1][MeterOverlap<double>::AVG_AVG], squared(stats[1][MeterOverlap<double>::AVG_ERR]));
    double blockCov = blockCor[0][1] * refAverage[threadNum].getStdDev() * refOverlapAverage.getStdDev();
    Uncertain<double>::setCovariance(refAverage[threadNum], refOverlapAverage, blockCov);
    refNumSteps[threadNum] = refMeter->getNumSamples();
    stats = targetMeter->getStatistics();
    targetAverage[threadNum]=Uncertain<double>(stats[0][MeterOverlap<double>::AVG_AVG], squared(stats[0][MeterOverlap<double>::AVG_ERR]));
    Uncertain<double> targetOverlapAverage = Uncertain<double>(stats[1][MeterOverlap<double>::AVG_AVG], squared(stats[1][MeterOverlap<double>::AVG_ERR]))/refMeter->getAlpha()[0];
    blockCor = targetMeter->getBlockCorrelation();
    blockCov = blockCor[0][1] * targetAverage[threadNum].getStdDev() * targetOverlapAverage.getStdDev();
    Uncertain<double>::setCovariance(targetAverage[threadNum], targetOverlapAverage, blockCov);
    targetNumSteps[threadNum] = targetMeter->getNumSamples();
}

void
ResultsVirial::
reduce() {
    if (reduced) {
        return;
    }
    refAverageReduced = 0.0;
    overlapRatioAverageReduced = 0.0;
    targetAverageReduced = 0.0;
    refNumStepsReduced = 0;
    targetNumStepsReduced = 0;
    virialCoefficientReduced = 0;

    for (int threadNum = 0; threadNum < numThreads; threadNum++) {
        refAverageReduced += refAverage[threadNum]/(double)numThreads;
        refNumStepsReduced += refNumSteps[threadNum];
        targetAverageReduced += targetAverage[threadNum]/(double)numThreads;
        targetNumStepsReduced += targetNumSteps[threadNum];
        overlapRatioAverageReduced += overlapRatioAverage[threadNum]/(double)numThreads;
        virialCoefficientReduced += virialCoefficient[threadNum]/(double)numThreads;
    }

#ifdef USE_MPI
    int mpiSize = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);

    const int mpiBufferSize = 10;

    double sendbuf[mpiBufferSize];

    int offset = 0;

    sendbuf[offset++] = refAverageReduced.getMean();
    sendbuf[offset++] = refAverageReduced.getVariance();
    sendbuf[offset++] = targetAverageReduced.getMean();
    sendbuf[offset++] = targetAverageReduced.getVariance();
    sendbuf[offset++] = refNumStepsReduced;
    sendbuf[offset++] = targetNumStepsReduced;
    sendbuf[offset++] = overlapRatioAverageReduced.getMean();
    sendbuf[offset++] = overlapRatioAverageReduced.getVariance();
    sendbuf[offset++] = virialCoefficientReduced.getMean();
    sendbuf[offset++] = virialCoefficientReduced.getVariance();

    double recvbuf[mpiBufferSize];

    for (int i = 0; i < mpiBufferSize; i++) {
        recvbuf[i] = 0;
    }

    MPI_Allreduce(sendbuf, recvbuf, mpiBufferSize, MPI_DOUBLE,
		  MPI_SUM, MPI_COMM_WORLD);

    offset = 0;

    refAverageReduced = Uncertain<double>(recvbuf[offset], recvbuf[offset+1]) / (double)mpiSize;
    offset += 2;
    targetAverageReduced = Uncertain<double>(recvbuf[offset], recvbuf[offset+1]) / (double)mpiSize;
    offset += 2;
    refNumStepsReduced = recvbuf[offset++];
    targetNumStepsReduced = recvbuf[offset++];
    overlapRatioAverageReduced = Uncertain<double>(recvbuf[offset], recvbuf[offset+1]) / (double)mpiSize;
    offset += 2;
    virialCoefficientReduced = Uncertain<double>(recvbuf[offset], recvbuf[offset+1]) / (double)mpiSize;
#endif
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
getTargetAverageReduced() const {
    assert(reduced);
    return targetAverageReduced;
}

Uncertain<double>
ResultsVirial::
getOverlapRatioAverageReduced() const {
    assert(reduced);
    return overlapRatioAverageReduced;
}

double
ResultsVirial::
getRefIntegral() const {
    return refIntegral;
}

void
ResultsVirial::
putOverlapRatio(int threadNum, double ratio, double uncertainty){
    overlapRatioAverage[threadNum] = Uncertain<double>(ratio, uncertainty);
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
