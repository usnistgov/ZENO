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
    bool verboseAdjust, tunable;

protected:
    std::string moveName;
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


/// Constructs the class to perform a monte carlo trial.
///

template <class T,
        class RandomNumberGenerator>
MCMove<T, RandomNumberGenerator>::
MCMove(IntegratorMSMC<T, RandomNumberGenerator> & integratorMSMC, ClusterSum<T> * clusterSum) : integratorMSMC(integratorMSMC), clusterSum(clusterSum)
              {
    stepSize = 0;
    numTrials = 0;
    numAccepted = 0;
    chiSum = 0;
    lastAdjust = 0;
    adjustInterval = 100;
    adjustStep = 1.05;
    minAdjustStep = 1;
    verboseAdjust = false;
    tunable = true;
    maxStepSize = 0;
}

template <class T,
        class RandomNumberGenerator>
MCMove<T, RandomNumberGenerator>::
  ~MCMove() {
}

/// Constructs a sub class of MCMove to perform a monte carlo trial for translation.
///
template <class T, 
        class RandomNumberGenerator>
MCMoveTranslate<T, RandomNumberGenerator>::
MCMoveTranslate(IntegratorMSMC<T, RandomNumberGenerator> & integratorMSMC, ClusterSum<T> * clusterSum) : MCMove<T, RandomNumberGenerator>(integratorMSMC, clusterSum)
{
    MCMove<T, RandomNumberGenerator>::moveName = "translate";
    MCMove<T, RandomNumberGenerator>::stepSize = cbrt(integratorMSMC.getParticles()->at(0)->numSpheres())*integratorMSMC.getParticles()->at(0)->getModel()->getSpheres()->at(0).getRadius();
    MCMove<T, RandomNumberGenerator>::maxStepSize = 1000;
}

template <class T,
        class RandomNumberGenerator>
MCMoveTranslate<T, RandomNumberGenerator>::
~MCMoveTranslate() {
}

/// Perform a monte carlo trial for translation.
///
template <class T,
        class RandomNumberGenerator>
double MCMoveTranslate<T, RandomNumberGenerator>::
doTrial(double oldValue, bool & accepted){
    MCMove<T, RandomNumberGenerator>::numTrials++;
    if (MCMove<T, RandomNumberGenerator>::tunable && MCMove<T, RandomNumberGenerator>::numTrials >= MCMove<T, RandomNumberGenerator>::adjustInterval) {
        MCMove<T, RandomNumberGenerator>::adjustStepSize();
    }
    if(oldValue == 0){
        std::cerr << "Old Value is zero " << std::endl;
        exit(1);
    }
    std::vector<Vector3<T>> step(MCMove<T, RandomNumberGenerator>::integratorMSMC.getParticles()->size() - 1);
    for(unsigned int j = 0; j < (MCMove<T, RandomNumberGenerator>::integratorMSMC.getParticles()->size() - 1); ++j)
    {
        for(int i = 0; i < 3; ++i)
        {
            step[j].set(i, (MCMove<T, RandomNumberGenerator>::integratorMSMC.getRandomNumberGenerator()->getRandIn01() - 0.5) * MCMove<T, RandomNumberGenerator>::stepSize);
        }
        MCMove<T, RandomNumberGenerator>::integratorMSMC.getParticles()->at(j+1)->translateBy(step[j]);
    }
    double newValue = MCMove<T, RandomNumberGenerator>::clusterSum->value();
    double ratio = newValue / oldValue;
    ratio = std::abs(ratio);
    MCMove<T, RandomNumberGenerator>::chiSum += ratio;
    accepted = (ratio > 1) || (ratio > MCMove<T, RandomNumberGenerator>::integratorMSMC.getRandomNumberGenerator()->getRandIn01());
    if(!accepted)
    {
        for(unsigned int j = 0; j < (MCMove<T, RandomNumberGenerator>::integratorMSMC.getParticles()->size() - 1); ++j)
        {
            step[j] *= -1;
            MCMove<T, RandomNumberGenerator>::integratorMSMC.getParticles()->at(j+1)->translateBy(step[j]);
        }
        newValue = oldValue;
    }
    return newValue;
}

/// Constructs a sub class of MCMove to perform a monte carlo trial for rotation.
///
template <class T,
        class RandomNumberGenerator>
MCMoveRotate<T, RandomNumberGenerator>::
MCMoveRotate(IntegratorMSMC<T, RandomNumberGenerator> & integratorMSMC, ClusterSum<T> * clusterSum) : MCMove<T, RandomNumberGenerator>(integratorMSMC, clusterSum)
{
    MCMove<T, RandomNumberGenerator>::moveName = "rotate";
    MCMove<T, RandomNumberGenerator>::stepSize = M_PI/4;
    MCMove<T, RandomNumberGenerator>::maxStepSize = M_PI/2;
}

template <class T, class RandomNumberGenerator>
MCMoveRotate<T, RandomNumberGenerator>::
~MCMoveRotate() {
}

/// Perform a monte carlo trial for rotation.
///
template <class T, class RandomNumberGenerator>
double MCMoveRotate<T, RandomNumberGenerator>::
doTrial(double oldValue, bool & accepted){
    MCMove<T, RandomNumberGenerator>::numTrials++;
    if (MCMove<T, RandomNumberGenerator>::tunable && MCMove<T, RandomNumberGenerator>::numTrials >= MCMove<T, RandomNumberGenerator>::adjustInterval) {
        MCMove<T, RandomNumberGenerator>::adjustStepSize();
    }
    if(oldValue == 0){
        std::cerr << "Old Value is zero " << std::endl;
        exit(1);
    }
    std::vector<Vector3<T>> axis(MCMove<T, RandomNumberGenerator>::integratorMSMC.getParticles()->size());
    std::vector<double> angle(MCMove<T, RandomNumberGenerator>::integratorMSMC.getParticles()->size());
    for(unsigned int j = 0; j < (MCMove<T, RandomNumberGenerator>::integratorMSMC.getParticles()->size()); ++j)
    {
        angle[j] =  (MCMove<T, RandomNumberGenerator>::integratorMSMC.getRandomNumberGenerator()->getRandIn01() - 0.5) * MCMove<T, RandomNumberGenerator>::stepSize;
        MCMove<T, RandomNumberGenerator>::integratorMSMC.getRandomUtilities()->setRandomOnSphere(&axis[j]);
        MCMove<T, RandomNumberGenerator>::integratorMSMC.getParticles()->at(j)->rotateBy(axis[j], angle[j]);
    }
    double newValue = oldValue;
    if(MCMove<T, RandomNumberGenerator>::clusterSum != NULL) {
        newValue = MCMove<T, RandomNumberGenerator>::clusterSum->value();
    }
    double ratio = newValue / oldValue;
    ratio = std::abs(ratio);
    MCMove<T, RandomNumberGenerator>::chiSum += ratio;

    accepted = (ratio > 1) || (ratio > MCMove<T, RandomNumberGenerator>::integratorMSMC.getRandomNumberGenerator()->getRandIn01());
    if(!accepted)
    {
        for(unsigned int j = 0; j < (MCMove<T, RandomNumberGenerator>::integratorMSMC.getParticles()->size()); ++j)
        {
            MCMove<T, RandomNumberGenerator>::integratorMSMC.getParticles()->at(j)->rotateBy(axis[j], -angle[j]);
        }
        newValue = oldValue;
    }
    return newValue;
}

/// Constructs a sub class of MCMove to perform a monte carlo trial for chain move.
///
template <class T, class RandomNumberGenerator>
MCMoveChainVirial<T, RandomNumberGenerator>::
MCMoveChainVirial(IntegratorMSMC<T, RandomNumberGenerator> & integratorMSMC, ClusterSum<T> * clusterSum, double sigma): MCMove<T, RandomNumberGenerator>(integratorMSMC, clusterSum), sigma(sigma)
{
    MCMove<T, RandomNumberGenerator>::moveName = "chain";
}

template <class T, class RandomNumberGenerator>
MCMoveChainVirial<T, RandomNumberGenerator>::
~MCMoveChainVirial() {
}

/// Perform a monte carlo trial for chain move.
///
template <class T, class RandomNumberGenerator>
double MCMoveChainVirial<T, RandomNumberGenerator>::
doTrial(double oldValue, bool & accepted) {
    MCMove<T, RandomNumberGenerator>::numTrials++;
    const Vector3<T> rPrev = MCMove<T, RandomNumberGenerator>::integratorMSMC.getParticles()->at(0)->getCenter();
    Vector3<T> sPrev = rPrev;
    accepted = true;
    for(unsigned int j = 1; j < MCMove<T, RandomNumberGenerator>::integratorMSMC.getParticles()->size(); ++j)
    {
        Vector3<T> r = MCMove<T, RandomNumberGenerator>::integratorMSMC.getParticles()->at(j)->getCenter();
        MCMove<T, RandomNumberGenerator>::integratorMSMC.getRandomUtilities()->setRandomInSphere(&r);
        Vector3<T> s = r*sigma + sPrev;
        MCMove<T, RandomNumberGenerator>::integratorMSMC.getParticles()->at(j)->setCenter(s);
         sPrev = s;
    }
    MCMove<T, RandomNumberGenerator>::chiSum += 1;
    return MCMove<T, RandomNumberGenerator>::clusterSum->value();
}

/// Adjusts step size.
///
template <class T, class RandomNumberGenerator>
void MCMove<T, RandomNumberGenerator>::
adjustStepSize(){
    double avg = chiSum/numTrials;
    if (avg > 0.5) {
        if (stepSize < maxStepSize) {
            if (lastAdjust < 0) {
                // back and forth
                adjustInterval *= 2;
                adjustStep = std::sqrt(adjustStep);
            }
            else if (lastAdjust == 5) {
                // sixth consecutive increase; increase adjustment step
                adjustStep *= adjustStep;
                if (adjustStep > 2) {
                    adjustStep = 2;
                }
                lastAdjust = 3;
            }
            stepSize *= adjustStep;
            stepSize = std::min(stepSize, maxStepSize);
            if (verboseAdjust) {
                printf("%s move increasing step size: %e (<chi> = %e)\n", MCMove<T, RandomNumberGenerator>::moveName.c_str(), stepSize, avg);
            }
            if (lastAdjust < 1) lastAdjust = 1;
            else lastAdjust++;
        }
    }
    else {
        if (lastAdjust > 0) {
            // back and forth
            adjustInterval *= 2;
            adjustStep = std::sqrt(adjustStep);
        }
        else if (lastAdjust == -5) {
            // sixth consecutive increase; increase adjustment step
            adjustStep *= adjustStep;
            if (adjustStep > 2) {
                adjustStep = 2;
            }
            lastAdjust = -3;
        }
        stepSize /= adjustStep;
        if (verboseAdjust) {
            printf("%s move decreasing step size: %e (<chi> = %e)\n", MCMove<T, RandomNumberGenerator>::moveName.c_str(), stepSize, avg);
        }
        if (lastAdjust > -1) lastAdjust = -1;
        else lastAdjust--;
    }
    numTrials = numAccepted = 0;
    chiSum = 0;
}

/// Returns step size.
///
template <class T, class RandomNumberGenerator>
double MCMove<T, RandomNumberGenerator>::
getStepSize() {
    return stepSize;
}

/// Sets step size.
///
template <class T, class RandomNumberGenerator>
void MCMove<T, RandomNumberGenerator>::
setStepSize(double sS) {
    stepSize = sS;
}

}
#endif
