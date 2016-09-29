// ================================================================
// 
// Disclaimer:  IMPORTANT:  This software was developed at the
// National Institute of Standards and Technology by employees of the
// Federal Government in the course of their official duties.
// Pursuant to title 17 Section 105 of the United States Code this
// software is not subject to copyright protection and is in the
// public domain.  This is an experimental system.  NIST assumes no
// responsibility whatsoever for its use by other parties, and makes
// no guarantees, expressed or implied, about its quality,
// reliability, or any other characteristic.  We would appreciate
// acknowledgement if the software is used.  This software can be
// redistributed and/or modified freely provided that any derivative
// works bear some notice that they are derived from it, and any
// modified versions bear some notice that they have been modified.
// 
// ================================================================

// ================================================================
// 
// Authors: Derek Juba <derek.juba@nist.gov>
// Date:    Fri Feb 13 13:31:22 2015 EDT
//
// Time-stamp: <2016-09-20 17:22:40 dcj>
//
// ================================================================

#ifndef SAMPLER_INTERIOR_H
#define SAMPLER_INTERIOR_H

#include <vector>

#include "../Geometry/Sphere.h"
#include "../Geometry/Vector3.h"

/// Samples random points inside a bounding sphere and determines whether they
/// hit an object, allowing for a given relative error in distance.
///
template <class T, 
  class RandomNumberGenerator, 
  class InsideOutsideTester,
  class RandomBallPointGenerator>
class SamplerInterior {
 public:
  SamplerInterior(RandomNumberGenerator * randomNumberGenerator, 
		  Sphere<T> const & boundingSphere, 
		  InsideOutsideTester const & insideOutsideTester,
		  T fracErrorBound);

  ~SamplerInterior();

  void sample(bool * hitObject,
	      Vector3<T> * hitPoint);

 private:
  RandomNumberGenerator * randomNumberGenerator;
  Sphere<T> const * boundingSphere;
  InsideOutsideTester const * insideOutsideTester;
  T fracErrorBound;
};

template <class T, 
  class RandomNumberGenerator, 
  class InsideOutsideTester,
  class RandomBallPointGenerator>
SamplerInterior<T, 
               RandomNumberGenerator, 
               InsideOutsideTester, 
               RandomBallPointGenerator>::
  SamplerInterior(RandomNumberGenerator * randomNumberGenerator, 
		 Sphere<T> const & boundingSphere, 
		 InsideOutsideTester const & insideOutsideTester,
		 T fracErrorBound) :
  randomNumberGenerator(randomNumberGenerator), 
  boundingSphere(&boundingSphere),
  insideOutsideTester(&insideOutsideTester), 
  fracErrorBound(fracErrorBound) {

}

template <class T, 
  class RandomNumberGenerator, 
  class InsideOutsideTester,
  class RandomBallPointGenerator>
SamplerInterior<T, 
               RandomNumberGenerator, 
               InsideOutsideTester, 
               RandomBallPointGenerator>::
  ~SamplerInterior() {

}

/// Compute a random point and determine whether it hits the object.
///
template <class T, 
  class RandomNumberGenerator, 
  class InsideOutsideTester,
  class RandomBallPointGenerator>
void 
SamplerInterior<T, 
               RandomNumberGenerator, 
               InsideOutsideTester, 
               RandomBallPointGenerator>::
  sample(bool * hitObject,
	 Vector3<T> * hitPoint) {

  Vector3<T> position = 
    RandomBallPointGenerator::generate(randomNumberGenerator, *boundingSphere);

  *hitObject = 
    insideOutsideTester->isInside(position,
				  fracErrorBound);

  *hitPoint = position;
}

#endif

// ================================================================

// Local Variables:
// time-stamp-line-limit: 30
// End:
