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
    virtual std::vector<double> doTrial(std::vector<double> oldValue, bool & accepted) = 0;
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

    std::vector<double> doTrial(std::vector<double> oldValue, bool & accepted);
};

/// Sub class of MCMove to perform a monte carlo trial for rotation.
///
template <class T,
        class RandomNumberGenerator>
class MCMoveRotate : public MCMove<T, RandomNumberGenerator> {
public:
    MCMoveRotate(IntegratorMSMC<T, RandomNumberGenerator> & integratorMSMC, ClusterSum<T> * clusterSum);
    ~MCMoveRotate();

    std::vector<double> doTrial(std::vector<double> oldValue, bool & accepted);
};

/// Sub class of MCMove to perform a monte carlo trial for chain move.
///
template <class T,
        class RandomNumberGenerator>
class MCMoveChainVirial : public MCMove<T, RandomNumberGenerator> {
public:
    MCMoveChainVirial(IntegratorMSMC<T, RandomNumberGenerator> & integratorMSMC, ClusterSum<T> * clusterSum, double sigma);
    ~MCMoveChainVirial();

    std::vector<double> doTrial(std::vector<double> oldValue, bool & accepted);
protected:
    double sigma;
};

struct StretchParameters {
  int sphere1;
  int sphere2;
  double step;
};

/// Sub class of MCMove to perform a monte carlo trial for bond stretch.
///
template <class T,
          class RandomNumberGenerator>
class MCMoveBondStretch : public MCMove<T, RandomNumberGenerator> {
public:
    MCMoveBondStretch(IntegratorMSMC<T, RandomNumberGenerator> & integratorMSMC, ClusterSum<T> * clusterSum, const Potential<T> & potential, double temperature);
    ~MCMoveBondStretch();

    std::vector<double> doTrial(std::vector<double> oldValue, bool & accepted);
protected:
    const Potential<T> & potential;
    double temperature;

    void transform(Vector3<T> dr, Particle<T> * p, int ip, Vector3<T> & shift);
    void transformBondedAtoms(Vector3<T> dr, Particle<T> * p, const std::vector<std::vector<BondedPartner>> * bondedPartners, int start, Vector3<T> & shift, std::vector<int> & modified);
};

struct AngleParameters {
  int sphere1;
  int sphere2;
  double step;
  Vector3<double> axis;
};

/// Sub class of MCMove to perform a monte carlo trial for bond angle bending.
///
template <class T,
          class RandomNumberGenerator>
class MCMoveBondAngle : public MCMove<T, RandomNumberGenerator> {
public:
    MCMoveBondAngle(IntegratorMSMC<T, RandomNumberGenerator> & integratorMSMC, ClusterSum<T> * clusterSum, const Potential<T> & potential, double temperature);
    ~MCMoveBondAngle();

    std::vector<double> doTrial(std::vector<double> oldValue, bool & accepted);
protected:
    const Potential<T> & potential;
    double temperature;

    void transform(Matrix3x3<T> & rotationTensor, Particle<T> * p, int ip, int b, Vector3<T> & shift);
    void transformBondedAtoms(Matrix3x3<T> & rotationTensor, Particle<T> * p, const std::vector<std::vector<BondedPartner>> * bondedPartners, int start, int b, Vector3<T> & shift, std::vector<int> & modified);
};

struct TorsionParameters {
  int sphere1;
  int sphere2;
  double step;
};

/// Sub class of MCMove to perform a monte carlo trial for bond angle bending.
///
template <class T,
          class RandomNumberGenerator>
class MCMoveBondTorsion : public MCMove<T, RandomNumberGenerator> {
public:
    MCMoveBondTorsion(IntegratorMSMC<T, RandomNumberGenerator> & integratorMSMC, ClusterSum<T> * clusterSum, const Potential<T> & potential, double temperature);
    ~MCMoveBondTorsion();

    std::vector<double> doTrial(std::vector<double> oldValue, bool & accepted);
protected:
    const Potential<T> & potential;
    double temperature;

    void transform(Matrix3x3<T> & rotationTensor, Particle<T> * p, int ip, int b, Vector3<T> & shift);
    void transformBondedAtoms(Matrix3x3<T> & rotationTensor, Particle<T> * p, const std::vector<std::vector<BondedPartner>> * bondedPartners, int start, int b, Vector3<T> & shift, std::vector<int> & modified);
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
std::vector<double> MCMoveTranslate<T, RandomNumberGenerator>::
doTrial(std::vector<double> oldValue, bool & accepted){
    MCMove<T, RandomNumberGenerator>::numTrials++;
    if (MCMove<T, RandomNumberGenerator>::tunable && MCMove<T, RandomNumberGenerator>::numTrials >= MCMove<T, RandomNumberGenerator>::adjustInterval) {
        MCMove<T, RandomNumberGenerator>::adjustStepSize();
    }
    if(oldValue[0] == 0){
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
    std::vector<double> newValue = MCMove<T, RandomNumberGenerator>::clusterSum->getValues();
    double ratio = newValue[0] / oldValue[0];
    ratio = std::abs(ratio);
    MCMove<T, RandomNumberGenerator>::chiSum += std::min(1.0, ratio);
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
std::vector<double> MCMoveRotate<T, RandomNumberGenerator>::
doTrial(std::vector<double> oldValue, bool & accepted){
    MCMove<T, RandomNumberGenerator>::numTrials++;
    if (MCMove<T, RandomNumberGenerator>::tunable && MCMove<T, RandomNumberGenerator>::numTrials >= MCMove<T, RandomNumberGenerator>::adjustInterval) {
        MCMove<T, RandomNumberGenerator>::adjustStepSize();
    }
    if(oldValue[0] == 0){
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
    std::vector<double> newValue = oldValue;
    if(MCMove<T, RandomNumberGenerator>::clusterSum != NULL) {
        newValue = MCMove<T, RandomNumberGenerator>::clusterSum->getValues();
    }
    double ratio = newValue[0] / oldValue[0];
    ratio = std::abs(ratio);
    MCMove<T, RandomNumberGenerator>::chiSum += std::min(1.0, ratio);

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
std::vector<double> MCMoveChainVirial<T, RandomNumberGenerator>::
doTrial(std::vector<double> oldValue, bool & accepted) {
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
    return MCMove<T, RandomNumberGenerator>::clusterSum->getValues();
}

/// Constructs a sub class of MCMove to perform a monte carlo trial for bond stretch move.
///
template <class T, class RandomNumberGenerator>
MCMoveBondStretch<T, RandomNumberGenerator>::
MCMoveBondStretch(IntegratorMSMC<T, RandomNumberGenerator> & integratorMSMC, ClusterSum<T> * clusterSum, const Potential<T> & potential, double temperature): MCMove<T, RandomNumberGenerator>(integratorMSMC, clusterSum), potential(potential), temperature(temperature)
{
    MCMove<T, RandomNumberGenerator>::moveName = "bond stretch";
    MCMove<T, RandomNumberGenerator>::stepSize = 0.1;
    MCMove<T, RandomNumberGenerator>::maxStepSize = 10;
}

template <class T, class RandomNumberGenerator>
MCMoveBondStretch<T, RandomNumberGenerator>::
~MCMoveBondStretch() {
}

/// Perform a monte carlo trial for chain move.
///
template <class T, class RandomNumberGenerator>
std::vector<double> MCMoveBondStretch<T, RandomNumberGenerator>::
doTrial(std::vector<double> oldValue, bool & accepted) {
    MCMove<T, RandomNumberGenerator>::numTrials++;
    if (MCMove<T, RandomNumberGenerator>::tunable && MCMove<T, RandomNumberGenerator>::numTrials >= MCMove<T, RandomNumberGenerator>::adjustInterval) {
        MCMove<T, RandomNumberGenerator>::adjustStepSize();
    }
    if(oldValue[0] == 0) {
        std::cerr << "Old Value is zero " << std::endl;
        exit(1);
    }
    double uOld = 0, uNew = 0;
    std::vector<StretchParameters> stepParameters;

    for(int j = 0; j < (int)MCMove<T, RandomNumberGenerator>::integratorMSMC.getParticles()->size(); j++) {
        Particle<T> * particle = MCMove<T, RandomNumberGenerator>::integratorMSMC.getParticles()->at(j);
        if (particle->numSpheres() < 2) continue;
        const std::vector<std::vector<BondedPartner>> * bondedPartners = potential.getBondedPartners(j);
        uOld += potential.energy1(j, particle);
        std::vector<int> modified;
        double s = (MCMove<T, RandomNumberGenerator>::integratorMSMC.getRandomNumberGenerator()->getRandIn01() - 0.5) * MCMove<T, RandomNumberGenerator>::stepSize;
        int b = 0;
        do {
            b = (int)(MCMove<T, RandomNumberGenerator>::integratorMSMC.getRandomNumberGenerator()->getRandIn01() * bondedPartners->size());
        }
        while (b == (int)bondedPartners->size() || (int)bondedPartners->at(b).size() < 2);
        modified.push_back(b);
        int a = 0;
        do {
            a = (int)(MCMove<T, RandomNumberGenerator>::integratorMSMC.getRandomNumberGenerator()->getRandIn01() * bondedPartners->at(b).size());
        }
        while (a == (int)bondedPartners->at(b).size());
        a = bondedPartners->at(b)[a].sphere2;
        modified.push_back(a);
        stepParameters.push_back({b, a, s});
        Vector3<T> dr = particle->getSpherePosition(b) - particle->getSpherePosition(a);
        dr *= s / dr.getMagnitude();
        Vector3<T> shift;
        transform(dr, particle, b, shift);
        transformBondedAtoms(dr, particle, bondedPartners, b, shift, modified);
        dr *= -1;
        transform(dr, particle, a, shift);
        transformBondedAtoms(dr, particle, bondedPartners, a, shift, modified);
        shift *= -1.0/particle->numSpheres();
        for (int k=0; k<particle->numSpheres(); k++) {
            particle->translateSphereBy(k, shift);
        }
        uNew += potential.energy1(j, particle);
    }
    std::vector<double> newValue = oldValue;
    if(MCMove<T, RandomNumberGenerator>::clusterSum != NULL) {
        newValue = MCMove<T, RandomNumberGenerator>::clusterSum->getValues();
    }
    double ratio = newValue[0] / oldValue[0] * std::exp(-(uNew-uOld)/temperature);
    ratio = std::abs(ratio);
    MCMove<T, RandomNumberGenerator>::chiSum += std::min(1.0, ratio);

    accepted = (ratio > 1) || (ratio > MCMove<T, RandomNumberGenerator>::integratorMSMC.getRandomNumberGenerator()->getRandIn01());
    if(!accepted)
    {
        for(unsigned int j = 0; j < (MCMove<T, RandomNumberGenerator>::integratorMSMC.getParticles()->size()); ++j)
        {
            Particle<T> * particle = MCMove<T, RandomNumberGenerator>::integratorMSMC.getParticles()->at(j);
            if (particle->numSpheres() < 2) continue;
            const std::vector<std::vector<BondedPartner>> * bondedPartners = potential.getBondedPartners(j);
            StretchParameters sp = stepParameters[j];
            int b = sp.sphere1, a = sp.sphere2;
            double s = -sp.step;
            std::vector<int> modified;
            modified.push_back(b);
            modified.push_back(a);
            Vector3<T> dr = particle->getSpherePosition(b) - particle->getSpherePosition(a);
            dr *= s / dr.getMagnitude();
            Vector3<T> shift;
            transform(dr, particle, b, shift);
            transformBondedAtoms(dr, particle, bondedPartners, b, shift, modified);
            dr *= -1;
            transform(dr, particle, a, shift);
            transformBondedAtoms(dr, particle, bondedPartners, a, shift, modified);
            shift *= -1.0/particle->numSpheres();
            for (int k=0; k<particle->numSpheres(); k++) {
              particle->translateSphereBy(k, shift);
            }
        }
        newValue = oldValue;
    }
    return newValue;
}

template <class T, class RandomNumberGenerator>
void MCMoveBondStretch<T, RandomNumberGenerator>::
transform(Vector3<T> dr, Particle<T> * p, int ip, Vector3<T> & shift) {
    p->translateSphereBy(ip, dr);
    shift += dr;
}

template <class T, class RandomNumberGenerator>
void MCMoveBondStretch<T, RandomNumberGenerator>::
transformBondedAtoms(Vector3<T> dr, Particle<T> * p, const std::vector<std::vector<BondedPartner>> * bondedPartners, int start, Vector3<T> & shift, std::vector<int> & modified) {
    for (int k=0; k<(int)bondedPartners->at(start).size(); k++) {
        bool rotated = false;
        int kk = bondedPartners->at(start)[k].sphere2;
        for (int ms : modified) {
            if (kk == ms) {
                rotated = true;
                break;
            }
        }
        if (!rotated) {
            transform(dr, p, kk, shift);
            modified.push_back(kk);
            transformBondedAtoms(dr, p, bondedPartners, kk, shift, modified);
        }
    }
}

/// Constructs a sub class of MCMove to perform a monte carlo trial for angle bending
///
template <class T, class RandomNumberGenerator>
MCMoveBondAngle<T, RandomNumberGenerator>::
MCMoveBondAngle(IntegratorMSMC<T, RandomNumberGenerator> & integratorMSMC, ClusterSum<T> * clusterSum, const Potential<T> & potential, double temperature): MCMove<T, RandomNumberGenerator>(integratorMSMC, clusterSum), potential(potential), temperature(temperature)
{
    MCMove<T, RandomNumberGenerator>::moveName = "bond angle";
    MCMove<T, RandomNumberGenerator>::stepSize = M_PI/8;
    MCMove<T, RandomNumberGenerator>::maxStepSize = M_PI/4;
}

template <class T, class RandomNumberGenerator>
MCMoveBondAngle<T, RandomNumberGenerator>::
~MCMoveBondAngle() {
}

/// Perform a monte carlo trial for angle bending move.
///
template <class T, class RandomNumberGenerator>
std::vector<double> MCMoveBondAngle<T, RandomNumberGenerator>::
doTrial(std::vector<double> oldValue, bool & accepted) {
  //std::cout << "angle doTrial" << std::endl;
    MCMove<T, RandomNumberGenerator>::numTrials++;
    if (MCMove<T, RandomNumberGenerator>::tunable && MCMove<T, RandomNumberGenerator>::numTrials >= MCMove<T, RandomNumberGenerator>::adjustInterval) {
        MCMove<T, RandomNumberGenerator>::adjustStepSize();
    }
    if(oldValue[0] == 0) {
        std::cerr << "Old Value is zero " << std::endl;
        exit(1);
    }
    double uOld = 0, uNew = 0;
    std::vector<AngleParameters> stepParameters;

    for(int j = 0; j < (int)MCMove<T, RandomNumberGenerator>::integratorMSMC.getParticles()->size(); j++) {
        Particle<T> * particle = MCMove<T, RandomNumberGenerator>::integratorMSMC.getParticles()->at(j);
        if (particle->numSpheres() < 3) continue;
        const std::vector<std::vector<BondedPartner>> * bondedPartners = potential.getBondedPartners(j);
        uOld += potential.energy1(j, particle);
        std::vector<int> modified;
        double s = (MCMove<T, RandomNumberGenerator>::integratorMSMC.getRandomNumberGenerator()->getRandIn01() - 0.5) * MCMove<T, RandomNumberGenerator>::stepSize;
        int b = 0;
        do {
            b = (int)(MCMove<T, RandomNumberGenerator>::integratorMSMC.getRandomNumberGenerator()->getRandIn01() * bondedPartners->size());
        }
        while (b == (int)bondedPartners->size() || (int)bondedPartners->at(b).size() < 2);
        modified.push_back(b);
        int a = 0;
        do {
            a = (int)(MCMove<T, RandomNumberGenerator>::integratorMSMC.getRandomNumberGenerator()->getRandIn01() * bondedPartners->at(b).size());
        }
        while (a == (int)bondedPartners->at(b).size());
        a = bondedPartners->at(b)[a].sphere2;
        modified.push_back(a);
        Vector3<double> axis;
        MCMove<T, RandomNumberGenerator>::integratorMSMC.getRandomUtilities()->setRandomOnSphere(&axis);
        Vector3<T> dr = particle->getSpherePosition(b) - particle->getSpherePosition(a);
        Vector3<T> projection = axis.dot(dr)/dr.getMagnitudeSqr() * dr;
        axis -= projection;
        axis.normalize();
        stepParameters.push_back({b, a, s, axis});
        Matrix3x3<T> rotation;
        rotation.setAxisAngle(axis, s);
        Vector3<T> shift;
        transform(rotation, particle, a, b, shift);
        transformBondedAtoms(rotation, particle, bondedPartners, a, b, shift, modified);
        shift *= -1.0/particle->numSpheres();
        for (int k=0; k<particle->numSpheres(); k++) {
            particle->translateSphereBy(k, shift);
        }
        uNew += potential.energy1(j, particle);
    }
    std::vector<double> newValue = oldValue;
    if(MCMove<T, RandomNumberGenerator>::clusterSum != NULL) {
        newValue = MCMove<T, RandomNumberGenerator>::clusterSum->getValues();
    }
    double ratio = newValue[0] / oldValue[0] * std::exp(-(uNew-uOld)/temperature);
    ratio = std::abs(ratio);
    MCMove<T, RandomNumberGenerator>::chiSum += std::min(1.0, ratio);

    accepted = (ratio > 1) || (ratio > MCMove<T, RandomNumberGenerator>::integratorMSMC.getRandomNumberGenerator()->getRandIn01());
    if(!accepted)
    {
        for(int j = 0; j < (int)(MCMove<T, RandomNumberGenerator>::integratorMSMC.getParticles()->size()); ++j)
        {
            Particle<T> * particle = MCMove<T, RandomNumberGenerator>::integratorMSMC.getParticles()->at(j);
            if (particle->numSpheres() < 3) continue;
            const std::vector<std::vector<BondedPartner>> * bondedPartners = potential.getBondedPartners(j);
            AngleParameters sp = stepParameters[j];
            int b = sp.sphere1, a = sp.sphere2;
            std::vector<int> modified;
            modified.push_back(b);
            modified.push_back(a);
            Matrix3x3<T> rotation;
            rotation.setAxisAngle(sp.axis, -sp.step);
            Vector3<T> shift;
            transform(rotation, particle, a, b, shift);
            transformBondedAtoms(rotation, particle, bondedPartners, a, b, shift, modified);
            shift *= -1.0/particle->numSpheres();
            for (int k=0; k<particle->numSpheres(); k++) {
                particle->translateSphereBy(k, shift);
            }
        }
        newValue = oldValue;
    }
    return newValue;
}

template <class T, class RandomNumberGenerator>
void MCMoveBondAngle<T, RandomNumberGenerator>::
transform(Matrix3x3<T> & rotationTensor, Particle<T> * p, int ip, int b, Vector3<T> & shift) {
    Vector3<T> dr = p->getSpherePosition(ip) - p->getSpherePosition(b);
    Vector3<T> dr0 = dr;
    rotationTensor.transform(dr);
    dr -= dr0;
    p->translateSphereBy(ip, dr);
    shift += dr;
}

template <class T, class RandomNumberGenerator>
void MCMoveBondAngle<T, RandomNumberGenerator>::
transformBondedAtoms(Matrix3x3<T> & rotationTensor, Particle<T> * p, const std::vector<std::vector<BondedPartner>> * bondedPartners, int start, int b, Vector3<T> & shift, std::vector<int> & modified) {
    for (int k=0; k<(int)bondedPartners->at(start).size(); k++) {
        bool rotated = false;
        int kk = bondedPartners->at(start)[k].sphere2;
        for (int ms : modified) {
            if (kk == ms) {
                rotated = true;
                break;
            }
        }
        if (!rotated) {
            transform(rotationTensor, p, kk, b, shift);
            modified.push_back(kk);
            transformBondedAtoms(rotationTensor, p, bondedPartners, kk, b, shift, modified);
        }
    }
}

/// Constructs a sub class of MCMove to perform a monte carlo trial for torsion angle rotation
///
template <class T, class RandomNumberGenerator>
MCMoveBondTorsion<T, RandomNumberGenerator>::
MCMoveBondTorsion(IntegratorMSMC<T, RandomNumberGenerator> & integratorMSMC, ClusterSum<T> * clusterSum, const Potential<T> & potential, double temperature): MCMove<T, RandomNumberGenerator>(integratorMSMC, clusterSum), potential(potential), temperature(temperature)
{
    MCMove<T, RandomNumberGenerator>::moveName = "bond angle";
    MCMove<T, RandomNumberGenerator>::stepSize = M_PI/8;
    MCMove<T, RandomNumberGenerator>::maxStepSize = M_PI/4;
}

template <class T, class RandomNumberGenerator>
MCMoveBondTorsion<T, RandomNumberGenerator>::
~MCMoveBondTorsion() {
}

/// Perform a monte carlo trial for angle bending move.
///
template <class T, class RandomNumberGenerator>
std::vector<double> MCMoveBondTorsion<T, RandomNumberGenerator>::
doTrial(std::vector<double> oldValue, bool & accepted) {
  //std::cout << "torsion doTrial" << std::endl;
    MCMove<T, RandomNumberGenerator>::numTrials++;
    if (MCMove<T, RandomNumberGenerator>::tunable && MCMove<T, RandomNumberGenerator>::numTrials >= MCMove<T, RandomNumberGenerator>::adjustInterval) {
        MCMove<T, RandomNumberGenerator>::adjustStepSize();
    }
    if(oldValue[0] == 0) {
        std::cerr << "Old Value is zero " << std::endl;
        exit(1);
    }
    double uOld = 0, uNew = 0;
    std::vector<TorsionParameters> stepParameters;

    for(int j = 0; j < (int)MCMove<T, RandomNumberGenerator>::integratorMSMC.getParticles()->size(); j++) {
        Particle<T> * particle = MCMove<T, RandomNumberGenerator>::integratorMSMC.getParticles()->at(j);
        if (particle->numSpheres() < 4) continue;
        const std::vector<std::vector<BondedPartner>> * bondedPartners = potential.getBondedPartners(j);
        uOld += potential.energy1(j, particle);
        std::vector<int> modified;
        double s = (MCMove<T, RandomNumberGenerator>::integratorMSMC.getRandomNumberGenerator()->getRandIn01() - 0.5) * MCMove<T, RandomNumberGenerator>::stepSize;
        int a = 0, b = 0;
        do {
            a = (int)(MCMove<T, RandomNumberGenerator>::integratorMSMC.getRandomNumberGenerator()->getRandIn01() * bondedPartners->size());
            if (a == (int)bondedPartners->size() || (int)bondedPartners->at(a).size() < 2) continue;
            b = (int)(MCMove<T, RandomNumberGenerator>::integratorMSMC.getRandomNumberGenerator()->getRandIn01() * bondedPartners->at(a).size());
            if (b == (int)bondedPartners->at(a).size()) continue;
            b = bondedPartners->at(a)[b].sphere2;
        } while ((int)bondedPartners->at(b).size() < 2);
        modified.push_back(a);
        modified.push_back(b);
        stepParameters.push_back({a, b, s});
        //std::cout << "torsion doTrial j " << j << " a " << a << " b " << b << std::endl;
        Vector3<double> axis = particle->getSpherePosition(b) - particle->getSpherePosition(a);
        axis.normalize();
        Matrix3x3<T> rotation;
        rotation.setAxisAngle(axis, s);
        Vector3<T> shift;
        // transform atoms bonded to a; axis passes through a
        transformBondedAtoms(rotation, particle, bondedPartners, a, a, shift, modified);
        // transform atoms bonded to b; axis passes through a
        transformBondedAtoms(rotation, particle, bondedPartners, b, a, shift, modified);
        shift *= -1.0/particle->numSpheres();
        for (int k=0; k<particle->numSpheres(); k++) {
            particle->translateSphereBy(k, shift);
        }
        uNew += potential.energy1(j, particle);
    }
    std::vector<double> newValue = oldValue;
    if(MCMove<T, RandomNumberGenerator>::clusterSum != NULL) {
        newValue = MCMove<T, RandomNumberGenerator>::clusterSum->getValues();
    }
    double ratio = newValue[0] / oldValue[0] * std::exp(-(uNew-uOld)/temperature);
    ratio = std::abs(ratio);
    MCMove<T, RandomNumberGenerator>::chiSum += std::min(1.0, ratio);

    accepted = (ratio > 1) || (ratio > MCMove<T, RandomNumberGenerator>::integratorMSMC.getRandomNumberGenerator()->getRandIn01());
  //std::cout << "torsion doTrial accepted? " << accepted << std::endl;
    if(!accepted)
    {
        for(int j = 0; j < (int)(MCMove<T, RandomNumberGenerator>::integratorMSMC.getParticles()->size()); ++j)
        {
            Particle<T> * particle = MCMove<T, RandomNumberGenerator>::integratorMSMC.getParticles()->at(j);
            if (particle->numSpheres() < 4) continue;
            const std::vector<std::vector<BondedPartner>> * bondedPartners = potential.getBondedPartners(j);
            TorsionParameters sp = stepParameters[j];
            int a = sp.sphere1, b = sp.sphere2;
            std::vector<int> modified;
            modified.push_back(a);
            modified.push_back(b);
            Vector3<double> axis = particle->getSpherePosition(b) - particle->getSpherePosition(a);
            axis.normalize();
            Matrix3x3<T> rotation;
            rotation.setAxisAngle(axis, -sp.step);
            Vector3<T> shift;
            transformBondedAtoms(rotation, particle, bondedPartners, a, a, shift, modified);
            transformBondedAtoms(rotation, particle, bondedPartners, a, b, shift, modified);
            shift *= -1.0/particle->numSpheres();
            for (int k=0; k<particle->numSpheres(); k++) {
                particle->translateSphereBy(k, shift);
            }
        }
        newValue = oldValue;
    }
  //std::cout << "torsion doTrial done" << std::endl;
    return newValue;
}

template <class T, class RandomNumberGenerator>
void MCMoveBondTorsion<T, RandomNumberGenerator>::
transform(Matrix3x3<T> & rotationTensor, Particle<T> * p, int ip, int b, Vector3<T> & shift) {
    Vector3<T> dr = p->getSpherePosition(ip) - p->getSpherePosition(b);
    Vector3<T> dr0 = dr;
    rotationTensor.transform(dr);
    dr -= dr0;
    p->translateSphereBy(ip, dr);
    shift += dr;
}

template <class T, class RandomNumberGenerator>
void MCMoveBondTorsion<T, RandomNumberGenerator>::
transformBondedAtoms(Matrix3x3<T> & rotationTensor, Particle<T> * p, const std::vector<std::vector<BondedPartner>> * bondedPartners, int start, int b, Vector3<T> & shift, std::vector<int> & modified) {
    for (int k=0; k<(int)bondedPartners->at(start).size(); k++) {
        bool rotated = false;
        int kk = bondedPartners->at(start)[k].sphere2;
        for (int ms : modified) {
            if (kk == ms) {
                rotated = true;
                break;
            }
        }
        if (!rotated) {
            transform(rotationTensor, p, kk, b, shift);
            modified.push_back(kk);
            transformBondedAtoms(rotationTensor, p, bondedPartners, kk, b, shift, modified);
        }
    }
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
