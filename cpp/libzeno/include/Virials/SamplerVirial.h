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

#ifndef SAMPLER_VIRIAL_H
#define SAMPLER_VIRIAL_H

#include <vector>

#include "../Timer.h"
#include "../Geometry/Sphere.h"
#include "../Geometry/MixedModel.h"
#include "OverlapTester.h"

namespace zeno {

/// Performs calculations to obtain virial coefficients.
template <class T,
  class RandomNumberGenerator>
class SamplerVirial {
 public:
  SamplerVirial(int threadNum,
                Timer const * totalTimer,
                RandomNumberGenerator * randomNumberGenerator,
                std::vector<Sphere<double> *> & boundingSpheres,
                std::vector<int> & numParticles,
                std::vector<MixedModel<T> *> & particles,
		        OverlapTester<T> const & overlapTester);

  ~SamplerVirial();

  void go(long long nSamples,
          double alpha,
          bool equilibrating,
          double refStepFrac);

 private:
  int threadNum;
  Timer const * totalTimer;
  RandomNumberGenerator * randomNumberGenerator;
  std::vector<Sphere<double> *> boundingSpheres;
  std::vector<int> & numParticles;
  std::vector<MixedModel<T> *> & particles;
  OverlapTester<T> const & overlapTester;
};

}

#endif

