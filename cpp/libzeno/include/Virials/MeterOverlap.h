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

#ifndef METER_OVERLAP_H
#define METER_OVERLAP_H

#include "ClusterSum.h"

///Collects data and statistics.
///
template <class T>
class ClusterSum;

template <class T>
class MeterOverlap {
 public:
    MeterOverlap(ClusterSum<T> & clusterSumPerturb, double alphaCenter, double alphaSpan, int numAlpha);
    ~MeterOverlap();
    void setAlpha(double alphaCenter, double alphaSpan);
    int getNumAlpha();
    const double * getAlpha(){return alpha;}
    void collectData(double primaryValue, bool accepted);
    double ** getStatistics() {
        if (blockCount == 0) {
            for (int i = 0; i < numData; ++i) {
                stats[i][AVG_CUR] = mostRecent[i];
                stats[i][1] = stats[i][2] = stats[i][3] = NAN;
            }
            return stats;
        }
        for (int i = 0; i < numData; ++i) {
            stats[i][AVG_CUR] = mostRecent[i];
            stats[i][AVG_AVG] = blockSum[i] / (blockSize * blockCount);
            if (blockCount == 1) {
                for (int i = 0; i < numData; ++i) {
                    stats[i][AVG_ERR] = stats[i][AVG_ACOR] = NAN;
                }
                continue;
            }
            stats[i][AVG_ERR] = blockSum2[i] / blockCount - stats[i][AVG_AVG] * stats[i][AVG_AVG];
            if (stats[i][AVG_ERR]<0) stats[i][AVG_ERR] = 0;
            if (stats[i][AVG_ERR] == 0) {
                stats[i][AVG_ACOR] = 0;
            }
            else {
                double bc;
                bc = (((2 * blockSum[i] / blockSize - firstBlockSum[i] - prevBlockSum[i]) * stats[i][AVG_AVG] - correlationSum[i]) / (1 - blockCount) + stats[i][AVG_AVG] * stats[i][AVG_AVG]) / stats[i][AVG_ERR];
                stats[i][AVG_ACOR] = (std::isnan(bc) || bc <= -1 || bc >= 1) ? 0 : bc;
            }
            stats[i][AVG_ERR] = sqrt(stats[i][AVG_ERR] / (blockCount - 1));
        }
        return stats;
    }

    /// Returns block covariance.
    ///

    double ** getBlockCovariance(){
        double totSamples = blockSize*blockCount;
        double totSq = totSamples*totSamples;
        for (int i = 0; i < numData; ++i) {
            for (int j = 0; j <= i; ++j) {
                blockCovariance[i][j] = blockCovariance[j][i] = blockCovSum[i][j] / blockCount - blockSum[i] * blockSum[j] / totSq;
            }
        }

        return blockCovariance;
    }

    /// Returns block correlation.
    ///
    double ** getBlockCorrelation(){
        getBlockCovariance();
        for (int i = 0; i < numData; ++i) {
            for (int j = 0; j < i; ++j) {
                double d = blockCovariance[i][i] * blockCovariance[j][j];
                double c = d <= 0 ? 0 : blockCovariance[i][j] / sqrt(d);
                blockCovariance[i][j] = blockCovariance[j][i] = c;
            }
        }
        for (int i = 0; i < numData; ++i) {
            blockCovariance[i][i] = 1;
        }
        return blockCovariance;
    }

    void setBlockSize(long long blockSize);
    long long getBlockSize() {return blockSize;}
    long long getBlockCount() {return blockCount;}
    long long getNumSamples() { return  blockCount*blockSize + blockSize - blockCountdown;}
    void setNumData(int newNumData);
    int getNumData();
    virtual void reset();
    double ** getRatioStatistics();
    double ** getRatioCovariance();
    double ** getRatioCorrelation();
    const static int AVG_CUR = 0, AVG_AVG = 1, AVG_ERR = 2, AVG_ACOR = 3;

    static double ratioErr(double nAvg, double nErr, double dAvg, double dErr, double cor) {
        if (nAvg == 0 && nErr == 0) return 0;
        double ratio = nAvg / dAvg;
        if (nAvg == 0) {
            return sqrt((nErr * nErr) / (dAvg * dAvg));
        }
        return sqrt((nErr * nErr / (nAvg * nAvg) + dErr * dErr / (dAvg * dAvg) - 2 * cor * nErr * dErr / (nAvg * dAvg)) * ratio * ratio);
    }
    /**
 * Compute covariance of i/d and j/d
 *
 * vi: value of i
 * vj: value of j
 * vd: value of d
 * ei: error in the i numerator
 * ej: error in the j numerator
 * ed: error in the denominator d
 * cij: correlation between i and j
 * cid: correlation between i and d
 * cjd: correlation between j and d
 */
    static double ratioCov(double vi, double vj, double vd, double ei, double ej, double ed, double cij, double cid, double cjd) {
        double eid = ratioErr(vi, ei, vd, ed, cid);
        double ejd = ratioErr(vj, ej, vd, ed, cjd);
        if (eid == 0 || ejd == 0) return 0;
        return (vi / vd) * (vj / vd) * ((ei / vi) * (ed / vd) + (ei / vi) * (ej / vj) * cij - (ei / vi) * (ed / vd) * cid - (ej / vj) * (ed / vj) * cjd);
    }

    /**
     * Compute correlation of i/d and j/d
     *
     * vi: value of i
     * vj: value of j
     * vd: value of d
     * ei: error in the i numerator
     * ej: error in the j numerator
     * ed: error in the denominator d
     * cij: correlation between i and j
     * cid: correlation between i and d
     * cjd: correlation between j and d
     */
    static double ratioCor(double vi, double vj, double vd, double ei, double ej, double ed, double cij, double cid, double cjd) {
        double eid = ratioErr(vi, ei, vd, ed, cid);
        double ejd = ratioErr(vj, ej, vd, ed, cjd);
        if (eid == 0 || ejd == 0) return 0;
        return ((vi / vd) / eid) * ((vj / vd) / ejd) * ((ed / vd) * (ed / vd) + (ei / vi) * (ej / vj) * cij - (ei / vi) * (ed / vd) * cid - (ej / vj) * (ed / vd) * cjd);
    }
private:
    ClusterSum<T> & clusterSumPerturb;
    double * data;
    double * alpha;
    const int numAlpha;
    int numData;
    long long blockSize, blockCount;
    long long blockCountdown;
    double * mostRecent;
    double * currentBlockSum, * blockSum, * blockSum2, * correlationSum;
    double * prevBlockSum, * firstBlockSum;
    double ** stats;
    double ** blockSums;
    double ** blockCovariance;
    double ** blockCovSum;
    double ** ratioStats;
    double ** ratioCovariance;
    double perturbValue;
};
#endif

