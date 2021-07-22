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

#include "Virials/MeterOverlap.h"
#include "Virials/alloc2D.h"

/// Constructs a class to collect data and statistics.
///

template <class T>
MeterOverlap<T>::
MeterOverlap(ClusterSum<T> & clusterSumPerturb,
             double alphaCenter, double alphaSpan, int numAlpha) :
             clusterSumPerturb(clusterSumPerturb),
             data(NULL), numAlpha(numAlpha), mostRecent(NULL), currentBlockSum(NULL), blockSum(NULL),
             blockSum2(NULL), correlationSum(NULL),prevBlockSum(NULL), firstBlockSum(NULL),
             stats(NULL), blockSums(NULL), blockCovariance(NULL), blockCovSum(NULL),
             ratioStats(NULL), ratioCovariance(NULL), perturbValue(-1){
    alpha = new double[numAlpha];
    setAlpha(alphaCenter, alphaSpan);
    setBlockSize(1000);
}

template <class T>
MeterOverlap<T>::
  ~MeterOverlap() {
    delete[] alpha;
    free(data);
    free(mostRecent);
    free(currentBlockSum);
    free(blockSum);
    free(blockSum2);
    free(correlationSum);
    free(prevBlockSum);
    free(firstBlockSum);
    free2D((void **)stats);
    free2D((void **)blockCovSum);
    free2D((void **)blockCovariance);
    free2D((void **)ratioStats);
    free2D((void **)ratioCovariance);
}

/// Sets alpha.
///
template <class T>
void
MeterOverlap<T>::
setAlpha(double alphaCenter, double alphaSpan) {
    if (numAlpha > 1) {
        numData = numAlpha;
        if (alphaSpan == 0) {
            std::cerr << "If # of alpha > 1, then alpha span can't be 0" << std::endl;
            exit(1);
        }
    }
    else {
        numData = 2;
        if (alphaSpan != 0) {
            std::cerr << "If # of alpha is 1, then alpha span must be 0" << std::endl;
            exit(1);
        }
    }

    if (numAlpha == 1) {
        alpha[0] = alphaCenter;
        reset();
        return;
    }
    for (int i = 0; i < numAlpha; ++i) {
        alpha[i] = alphaCenter * exp((i - (numAlpha - 1.0) / 2) / (numAlpha - 1.0) * alphaSpan);
    }
    reset();
}

/// Returns numAlpha.
///
template <class T>
int
MeterOverlap<T>::
getNumAlpha() {
 return numAlpha;
}

/// Collects data.
///
template <class T>
void
MeterOverlap<T>::
collectData(double primaryValue, bool accepted) {
    double pi = fabs(primaryValue);
    if (pi == 0 || pi == std::numeric_limits<double>::infinity() || std::isnan(pi)) {
        std::cerr << "pi is" << pi << std::endl;
        exit(1);
    }
    if(accepted || perturbValue == -1){
        perturbValue = fabs(clusterSumPerturb.value());
    }
    if (numAlpha == 1) {
        data[0] = primaryValue / pi;
        data[1] = perturbValue / (perturbValue + alpha[0] * pi);
    } else {
        for (int i = 0; i < numAlpha; ++i) {
            // gamma_OS = pi1 pi0 / (pi1 + alpha pi0)
            // 0: gamma_OS/pi0 = pi1 / (pi1 + alpha pi0)
            // 1: gamma_OS/pi1 = pi0 / (pi1 + alpha pi0)
            //                 = (1/alpha) pi0 / (pi1 + (1/alpha) pi0)
            // for 1 case, we use negative alphaSpan (alpha => 1/alpha)
            //    and effectively compute: alpha gammaOS/pi1
            // <0>/<1> = (1/alpha) <gammaOS/pi0>0 / <gammaOS/pi1>1
            //         ~= 1 (when alpha is optimal)
            data[i] = perturbValue / (perturbValue + alpha[i] * pi);
        }
    }
    for (int i = 0; i < numData; ++i) {
        mostRecent[i] = data[i];
        currentBlockSum[i] += data[i];
    }
    if (--blockCountdown == 0) {
        double blockSizeSq = ((double) blockSize) * ((double) blockSize);
        for (int i = 0; i < numData; ++i) {
            for (int j = 0; j <= i; ++j) {
                double ijx = currentBlockSum[i] * currentBlockSum[j] / blockSizeSq;
                blockCovSum[i][j] += ijx;
            }
        }
        for (int i = 0; i < numData; ++i) {
            blockSum[i] += currentBlockSum[i];
            currentBlockSum[i] /= blockSize;
            if (blockCount > 0) correlationSum[i] += prevBlockSum[i] * currentBlockSum[i];
            else firstBlockSum[i] = currentBlockSum[i];
            prevBlockSum[i] = currentBlockSum[i];
            blockSum2[i] += currentBlockSum[i] * currentBlockSum[i];
            currentBlockSum[i] = 0;
        }
        blockCount++;
        blockCountdown = blockSize;
    }
}

/// Sets block size.
///
template <class T>
void
MeterOverlap<T>::
setBlockSize(long long bs) {
    blockSize = bs;
    reset();
}

/// Resets the averages and statistics.
///
template <class T>
void
MeterOverlap<T>::
reset() {
    blockCount = 0;
    blockCountdown = blockSize;
    if (numData == 0) {
        // we can't allocate 0-size arrays, so just leave them as nullptr
        // at some point numData will be positive, reset will be called again
        return;
    }
    // realloc our arrays so that we can adjust if n changes
    data = (double *) realloc(data, numData * sizeof(double));
    mostRecent = (double *) realloc(mostRecent, numData * sizeof(double));
    currentBlockSum = (double *) realloc(currentBlockSum, numData * sizeof(double));
    blockSum = (double *) realloc(blockSum, numData * sizeof(double));
    blockSum2 = (double *) realloc(blockSum2, numData * sizeof(double));
    correlationSum = (double *) realloc(correlationSum, numData * sizeof(double));
    for (int i = 0; i < numData; ++i) {
        currentBlockSum[i] = blockSum[i] = blockSum2[i] = correlationSum[i] = 0;
    }
    prevBlockSum = (double *) realloc(prevBlockSum, numData * sizeof(double));
    firstBlockSum = (double *) realloc(firstBlockSum, numData * sizeof(double));
    stats = (double **) realloc2D((void **) stats, numData, 4, sizeof(double));
    blockCovSum = (double **) realloc2D((void **) blockCovSum, numData, numData, sizeof(double));
    for (int i = 0; i < numData; ++i) {
        for (int j = 0; j <= i; ++j) {
            blockCovSum[i][j] = 0;
        }
    }
    blockCovariance = (double **) realloc2D((void **) blockCovariance, numData, numData, sizeof(double));
    ratioStats = (double ** ) realloc2D ((void ** ) ratioStats, numData, 3, sizeof(double));
    ratioCovariance = (double ** ) realloc2D ((void ** ) ratioCovariance, numData, numData, sizeof(double));
}

/// Returns ratio statistics.
///
template <class T>
double **
MeterOverlap<T>::
getRatioStatistics() {
    if (blockCount == 0) {
        for (int i = 0; i < numData; ++i) {
            ratioStats[i][AVG_CUR] = NAN;
            ratioStats[i][AVG_AVG] = ratioStats[i][AVG_ERR] = NAN;
        }
        return ratioStats;
    }
    getStatistics();
    getBlockCovariance();
    for (int i = 0; i < numData; ++i) {
        ratioStats[i][AVG_CUR] = mostRecent[i] / mostRecent[numData - 1];
        ratioStats[i][AVG_AVG] = blockSum[i] / blockSum[numData - 1];
        if (blockCount == 1) {
            for (int i = 0; i < numData; ++i) {
                ratioStats[i][AVG_ERR] = NAN;
            }
            continue;
        }
        double d = blockCovariance[i][i] * blockCovariance[numData - 1][numData - 1];
        double icor = d <= 0 ? 0 : blockCovariance[i][numData - 1] / sqrt(d);
        ratioStats[i][AVG_ERR] = ratioErr(stats[i][AVG_AVG], stats[i][AVG_ERR], stats[numData - 1][AVG_AVG], stats[numData - 1][AVG_ERR], icor);
    }
    return ratioStats;
}

/// Returns ratio covariance.
///
template <class T>
double **
MeterOverlap<T>::
getRatioCovariance() {
    if (blockCount<2) {
        for (int i = 0; i < numData; ++i) {
            for (int j = 0; j < numData; ++j) {
                ratioCovariance[i][j] = NAN;
            }
        }
        return ratioCovariance;
    }
    getStatistics();
    getBlockCovariance();
    double vd = stats[numData - 1][AVG_AVG];
    double ed = stats[numData - 1][AVG_ERR];
    for (int i = 0; i < numData; ++i) {
        double vi = stats[i][AVG_AVG];
        double ei = stats[i][AVG_ERR];
        double x = blockCovariance[i][i] * blockCovariance[numData - 1][numData - 1];
        if (x <= 0) {
            for (int j = 0; j <= i; ++j) ratioCovariance[j][i] = ratioCovariance[i][j] = 0;
            continue;
        }
        double cid = blockCovariance[i][numData - 1] / sqrt(x);
        for (int j = 0; j <= i; ++j) {
            double vj = stats[j][AVG_AVG];
            double ej = stats[j][AVG_ERR];
            if (blockCovariance[j][j] <= 0) {
                ratioCovariance[j][i] = ratioCovariance[i][j] = 0;
                continue;
            }
            double cjd = blockCovariance[j][numData - 1] / sqrt(blockCovariance[j][j] * blockCovariance[numData - 1][numData - 1]);
            double cij = blockCovariance[i][j] / sqrt(blockCovariance[i][i]*blockCovariance[j][j]);
            ratioCovariance[j][i] = ratioCovariance[i][j] = ratioCov(vi, vj, vd, ei, ej, ed, cij, cid, cjd);
        }
    }
    return ratioCovariance;
}

/// Returns ratio correlation.
///
template <class T>
double **
MeterOverlap<T>::
getRatioCorrelation() {
    if (blockCount<2) {
        for (int i = 0; i < numData; ++i) {
            for (int j = 0; j < numData; ++j) {
                ratioCovariance[i][j] = NAN;
            }
        }
        return ratioCovariance;
    }
    getStatistics();
    getBlockCovariance();
    double vd = stats[numData - 1][AVG_AVG];
    double ed = stats[numData - 1][AVG_ERR];
    for (int i = 0; i < numData; ++i) {
        double vi = stats[i][AVG_AVG];
        double ei = stats[i][AVG_ERR];
        double x = blockCovariance[i][i] * blockCovariance[numData - 1][numData - 1];
        if (x <= 0) {
            for (int j = 0; j <= i; ++j) ratioCovariance[j][i] = ratioCovariance[i][j] = 0;
            continue;
        }
        double cid = blockCovariance[i][numData - 1] / sqrt(x);
        for (int j = 0; j <= i; ++j) {
            double vj = stats[j][AVG_AVG];
            double ej = stats[j][AVG_ERR];
            if (blockCovariance[j][j] <= 0) {
                ratioCovariance[j][i] = ratioCovariance[i][j] = 0;
                continue;
            }
            double cjd = blockCovariance[j][numData - 1] / sqrt(blockCovariance[j][j] * blockCovariance[numData - 1][numData - 1]);
            double cij = blockCovariance[i][j] / sqrt(blockCovariance[i][i]*blockCovariance[j][j]);
            ratioCovariance[j][i] = ratioCovariance[i][j] = ratioCor(vi, vj, vd, ei, ej, ed, cij, cid, cjd);
        }
    }
    return ratioCovariance;
}

template<class T>
int MeterOverlap<T>::getNumData() {
    return numData;
}
