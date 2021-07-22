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

#include "Virials/SamplerVirial.h"

///  Constructs the class to perform calculations to obtain virial coefficients.

using namespace zeno;

template <class T,
  class RandomNumberGenerator>
SamplerVirial<T,
               RandomNumberGenerator>::
  SamplerVirial(int threadNum,
                Timer const * totalTimer,
                RandomNumberGenerator * randomNumberGenerator,
                std::vector<Sphere<double> *> & boundingSpheres,
                std::vector<int> & numParticles,
                std::vector<MixedModel<T> *> & particles,
                OverlapTester<T> const & overlapTester) :
              threadNum(threadNum),
              totalTimer(totalTimer),
              randomNumberGenerator(randomNumberGenerator),
              boundingSpheres(boundingSpheres),
              numParticles(numParticles),
              particles(particles),
              overlapTester(overlapTester) {

}

template <class T,
  class RandomNumberGenerator>
SamplerVirial<T,
               RandomNumberGenerator>::
  ~SamplerVirial() {

}

/// Computes something.
///
template <class T,
  class RandomNumberGenerator>
void
SamplerVirial<T,
               RandomNumberGenerator>::
  go(long long nSamples,
     double alpha,
     bool equilibrating,
     double refStepFrac) {

    int numTargetBlocks = 0, numReferenceBlocks = 0;
    int nBlocks = 1000;

    if (nSamples < 100)
    {
        nBlocks = 1;
    }
    else if(nSamples < 1000)
    {
        nBlocks = 10;
    }
    else if(nSamples < 10000)
    {
        nBlocks = 100;
    }

   for(int step = 0; step < nBlocks; ++step)
   {
       bool runTarget = step*refStepFrac < numReferenceBlocks;
       for(long long subStep = 0; subStep < nSamples / nBlocks; ++subStep)
       {

       }
       if(runTarget) ++numTargetBlocks;
       else ++numReferenceBlocks;
   }
}

