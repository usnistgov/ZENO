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

#include "Virials/VirialAlpha.h"

// perhaps just take Clusters and make Meter and Average internally

template <class T,
        class RandomNumberGenerator>
VirialAlpha<T, RandomNumberGenerator>::
VirialAlpha(IntegratorMSMC<T, RandomNumberGenerator> & rIntegrator,
            IntegratorMSMC<T, RandomNumberGenerator> & tIntegrator,
            ClusterSum<T> & refClusterRef,
            ClusterSum<T> & refClusterTarget,
            ClusterSum<T> & targetClusterRef,
            ClusterSum<T> & targetClusterTarget) :
            stepCount(0), nextCheck(1000), refIntegrator(rIntegrator), targetIntegrator(tIntegrator),
            refMeter(MeterOverlap<T>(refClusterTarget, 1, 5, 10)),
            targetMeter(MeterOverlap<T>(targetClusterRef, 1, -5, 10)),
            newAlpha(0), newAlphaErr(0), alphaCor(0), alphaSpan(0), allDone(false),
            verbose(false) {
    int numAlpha = refMeter.getNumAlpha();
    const double *alpha = refMeter.getAlpha();
    alphaSpan = log(alpha[numAlpha-1]/alpha[0]);
    refIntegrator.setMeter(&refMeter);
    targetIntegrator.setMeter(&targetMeter);
}

template <class T,
        class RandomNumberGenerator>
VirialAlpha<T, RandomNumberGenerator>::
  ~VirialAlpha() {
}

template <class T,
        class RandomNumberGenerator>
void
VirialAlpha<T, RandomNumberGenerator>::
setVerbose(bool newVerbose) {
    verbose = newVerbose;
}

template <class T,
        class RandomNumberGenerator>
void
VirialAlpha<T, RandomNumberGenerator>::
setAlpha(double aCenter, double aSpan) {
    alphaSpan = aSpan;
    refMeter.setAlpha(aCenter, aSpan);
    refMeter.reset();
    targetMeter.setAlpha(1/aCenter, -aSpan);
    targetMeter.reset();
}

template <class T,
        class RandomNumberGenerator>
void
VirialAlpha<T, RandomNumberGenerator>::
analyze(double &jBest) {
    newAlpha = 0;
    newAlphaErr = 0;
    int numAlpha = refMeter.getNumAlpha();
    double lnRatio[numAlpha];
    const double *alpha = refMeter.getAlpha();
    double **refOverStats = refMeter.getStatistics();
    double **targetOverStats = targetMeter.getStatistics();
    for (int j=0; j<numAlpha; j++) {
        lnRatio[j] = log(refOverStats[j][MeterOverlap<T>::AVG_AVG]/targetOverStats[j][MeterOverlap<T>::AVG_AVG]);
        if (j>0 && lnRatio[j]*lnRatio[j-1] <= 0) {
            // linear interpolation on log scale
            double xj = lnRatio[j-1]/(lnRatio[j-1]-lnRatio[j]);
            jBest = j-1 + xj;
            newAlpha = exp(log(alpha[j-1]) + xj*(log(alpha[j]/alpha[j-1])));

            double ratio1 = exp(lnRatio[j-1]);
            double err1 = MeterOverlap<T>::ratioErr(refOverStats[j-1][MeterOverlap<T>::AVG_AVG], refOverStats[j-1][MeterOverlap<T>::AVG_ERR],
                                                 targetOverStats[j-1][MeterOverlap<T>::AVG_AVG], targetOverStats[j-1][MeterOverlap<T>::AVG_ERR], 0);
            double ratio2 = exp(lnRatio[j]);
            double err2 = MeterOverlap<T>::ratioErr(refOverStats[j][MeterOverlap<T>::AVG_AVG], refOverStats[j][MeterOverlap<T>::AVG_ERR],
                                                 targetOverStats[j][MeterOverlap<T>::AVG_AVG], targetOverStats[j][MeterOverlap<T>::AVG_ERR], 0);
            newAlphaErr = (err1/ratio1 > err2/ratio2 ? err1/ratio1 : err2/ratio2)*newAlpha;
            double ac1 = targetOverStats[j-1][MeterOverlap<T>::AVG_ACOR];
            double ac2 = targetOverStats[j][MeterOverlap<T>::AVG_ACOR];
            alphaCor = ac1 > ac2 ? ac1 : ac2;
            return;
        }
    }
    int jb = (fabs(lnRatio[0]) < fabs(lnRatio[numAlpha-1])) ? 0 : numAlpha-1;
    jBest = jb;
    newAlpha = alpha[jb];
    newAlphaErr = MeterOverlap<T>::ratioErr(refOverStats[jb][MeterOverlap<T>::AVG_AVG], refOverStats[jb][MeterOverlap<T>::AVG_ERR],
                                         targetOverStats[jb][MeterOverlap<T>::AVG_AVG], targetOverStats[jb][MeterOverlap<T>::AVG_ERR], 0);
}

template <class T,
        class RandomNumberGenerator>
void
VirialAlpha<T, RandomNumberGenerator>::
getNewAlpha(double &na, double &nae, double &ac) {
    na = newAlpha;
    nae = newAlphaErr;
    ac = alphaCor;
}

template <class T,
        class RandomNumberGenerator>
double *
VirialAlpha<T, RandomNumberGenerator>::
getAlphaStatistics() {
    alphaStats[0] = newAlpha;
    alphaStats[1] = newAlphaErr;
    alphaStats[2] = alphaCor;
    alphaStats[3] = alphaSpan;
    return alphaStats;
}

template <class T,
        class RandomNumberGenerator>
void
VirialAlpha<T, RandomNumberGenerator>::
run() {
    while (!allDone) {
        runSteps(1000);
    }
}

template <class T,
        class RandomNumberGenerator>
bool
VirialAlpha<T, RandomNumberGenerator>::
getAllDone() {
    return allDone;
}

template <class T,
        class RandomNumberGenerator>
void
VirialAlpha<T, RandomNumberGenerator>::
runSteps(int numSteps) {
    refIntegrator.doStep(numSteps);
    targetIntegrator.doStep(numSteps);
    stepCount += numSteps;
    if (stepCount >= nextCheck) {
        double jBestAlpha;
        analyze(jBestAlpha);
        if (verbose) printf("alpha  avg: %22.15e   err: %12.5e   cor: % 6.4f\n", newAlpha, newAlphaErr, alphaCor);
        int numAlpha = refMeter.getNumAlpha();
        double nextCheckFac = 1.4;
        if (jBestAlpha<numAlpha*0.1 || jBestAlpha>(numAlpha-1)*0.9) alphaSpan *= 2;
        else if (alphaCor < 0.3 && alphaSpan > 0.5 && jBestAlpha>numAlpha*0.2 && jBestAlpha<(numAlpha-1)*0.8) alphaSpan *= 0.25;
        else if (alphaCor < 0.6 && alphaSpan > 0.5 && jBestAlpha>numAlpha*0.2 && jBestAlpha<(numAlpha-1)*0.8) alphaSpan *= 0.6;
        else if (alphaCor > 0.2) nextCheckFac *= 2;
        else if (alphaCor < 0.1 && newAlphaErr/newAlpha < 0.02) allDone = true;
        setAlpha(newAlpha, alphaSpan);
        nextCheck *= nextCheckFac;
        nextCheck = stepCount + nextCheck;
    }
}
