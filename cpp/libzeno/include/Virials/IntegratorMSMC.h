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

#ifndef INTEGRATOR_MSMC_H
#define INTEGRATOR_MSMC_H

#include <vector>

#include "../Timer.h"
#include "../Geometry/Sphere.h"
#include "../Geometry/MixedModel.h"
#include "OverlapTester.h"
#include "MCMove.h"
#include "MeterOverlap.h"
#include "Particle.h"
#include "RandomUtilities.h"

namespace zeno {

/// Performs a step. A step performed by the integrator consists of selecting a MCMove, performing the trial defined by MCMove and collecting data and statistics.
///

template <class T,
  class RandomNumberGenerator>
class IntegratorMSMC {
 public:
    IntegratorMSMC(int threadNum,
                   Timer const * totalTimer,
                   RandomNumberGenerator * randomNumberGenerator,
                std::vector<Sphere<double> const *> & boundingSpheres,
                std::vector<int> & numParticles,
                std::vector<MixedModel<T> const *> & models);

  ~IntegratorMSMC();

  void doStep(long long numSteps);
  std::vector<Particle<T> *> * getParticles();
  RandomNumberGenerator * getRandomNumberGenerator();
  RandomUtilities<T, RandomNumberGenerator> * getRandomUtilities();
  void addMove(MCMove<T, RandomNumberGenerator> * mcMove, double moveProb);
  void setMeter(MeterOverlap<T> * meter);
  void setCurrentValue(double currentValue);
private:
  int threadNum;
  Timer const * totalTimer;
  RandomNumberGenerator * randomNumberGenerator;
  std::vector<Sphere<double> const *> & boundingSpheres;
  std::vector<int> & numParticles;
  std::vector<Particle<T> *> particles;
  std::vector<MCMove<T, RandomNumberGenerator> *> mcMoves;
  std::vector<double> moveProbs;
  RandomUtilities<T, RandomNumberGenerator> randomUtilities;
  MeterOverlap<T> * meterOverlap;
  double currentValue;
};

}
#endif

