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

#include "Virials/IntegratorMSMC.h"
#include "Virials/MCMove.h"

///  Constructs the class to perform a step.
///

using namespace zeno;

template <class T,
  class RandomNumberGenerator>
IntegratorMSMC<T,
               RandomNumberGenerator>::
IntegratorMSMC(int threadNum,
                Timer const * totalTimer,
                RandomNumberGenerator * randomNumberGenerator,
                std::vector<Sphere<double> const *> & boundingSpheres,
                std::vector<int> & numParticles,
                std::vector<MixedModel<T> const *> & models) :
              threadNum(threadNum),
              totalTimer(totalTimer),
              randomNumberGenerator(randomNumberGenerator),
              boundingSpheres(boundingSpheres),
              numParticles(numParticles), randomUtilities(randomNumberGenerator), currentValue(0){
    for(unsigned int i = 0; i < numParticles.size(); ++i)
    {
        for (int j =0; j < numParticles[i]; ++j)
        {
            particles.push_back(new Particle<T> (* models[i], * boundingSpheres[i]));
        }
    }

}

template <class T,
  class RandomNumberGenerator>
IntegratorMSMC<T,
               RandomNumberGenerator>::
  ~IntegratorMSMC() {

}

/// Carries out a step. A step performed by the integrator consists of selecting a MCMove, performing the trial defined by MCMove and collecting data and statistics.
///
template <class T,
  class RandomNumberGenerator>
void
IntegratorMSMC<T,
               RandomNumberGenerator>::
doStep(long long numSteps){
    double totalProb = 0.0;
    for(unsigned int i = 0; i < moveProbs.size(); ++i) totalProb += moveProbs[i];
    for(int j = 0; j < numSteps; ++j)
    {
        bool accepted = false;
        double cumProb = 0.0;
        double random = totalProb * randomNumberGenerator->getRandIn01();
        for(unsigned int i = 0; i < mcMoves.size(); ++i)
        {
            cumProb += moveProbs[i];
            if(random < cumProb)
            {
                currentValue = mcMoves[i]->doTrial(currentValue, accepted);
                break;
            }
        }
        meterOverlap->collectData(currentValue, accepted);
    }
}

/// Returns particles.
///
template <class T,
        class RandomNumberGenerator>
std::vector<Particle<T> *> *
IntegratorMSMC<T,RandomNumberGenerator>::
getParticles(){
    return &particles;
}

/// Returns random number generator.
///
template <class T,
        class RandomNumberGenerator>
RandomNumberGenerator *
IntegratorMSMC<T,RandomNumberGenerator>::
getRandomNumberGenerator(){
    return randomNumberGenerator;
}

/// Returns random utilities.
///
template <class T,
        class RandomNumberGenerator>
RandomUtilities<T, RandomNumberGenerator> *
IntegratorMSMC<T,RandomNumberGenerator>::
getRandomUtilities(){
    return & randomUtilities;
}


/// Adds a MC Move with Move Probability.
///
template <class T,
        class RandomNumberGenerator>
void
IntegratorMSMC<T, RandomNumberGenerator>::
addMove(MCMove<T, RandomNumberGenerator> * mcMove, double moveProb){
    mcMoves.push_back(mcMove);
    moveProbs.push_back(moveProb);
}

/// Adds a meter.
///
template <class T,
        class RandomNumberGenerator>
void
IntegratorMSMC<T, RandomNumberGenerator>::
setMeter(MeterOverlap<T> * meter) {
    meterOverlap = meter;
}

/// Adds current cluster value.
///
template <class T,
        class RandomNumberGenerator>
void
IntegratorMSMC<T, RandomNumberGenerator>::
setCurrentValue(double currentValue) {
    this->currentValue = currentValue;
}
