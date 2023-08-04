// ================================================================
//
// This software was developed by employees of the National Institute of
// Standards and Technology (NIST), an agency of the Federal Government.
// Pursuant to title 17 United States Code Section 105, works of NIST employees
// are not subject to copyright protection in the United States and are
// considered to be in the public domain. Permission to freely use, copy,
// modify, and distribute this software and its documentation without fee is
// hereby granted, provided that this notice and disclaimer of warranty appears
// in all copies.
//
// THE SOFTWARE IS PROVIDED 'AS IS' WITHOUT ANY WARRANTY OF ANY KIND, EITHER
// EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, ANY
// WARRANTY THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED
// WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND FREEDOM
// FROM INFRINGEMENT, AND ANY WARRANTY THAT THE DOCUMENTATION WILL CONFORM TO
// THE SOFTWARE, OR ANY WARRANTY THAT THE SOFTWARE WILL BE ERROR FREE. IN NO
// EVENT SHALL NIST BE LIABLE FOR ANY DAMAGES, INCLUDING, BUT NOT LIMITED TO,
// DIRECT, INDIRECT, SPECIAL OR CONSEQUENTIAL DAMAGES, ARISING OUT OF,
// RESULTING FROM, OR IN ANY WAY CONNECTED WITH THIS SOFTWARE, WHETHER OR NOT
// BASED UPON WARRANTY, CONTRACT, TORT, OR OTHERWISE, WHETHER OR NOT INJURY WAS
// SUSTAINED BY PERSONS OR PROPERTY OR OTHERWISE, AND WHETHER OR NOT LOSS WAS
// SUSTAINED FROM, OR AROSE OUT OF THE RESULTS OF, OR USE OF, THE SOFTWARE OR
// SERVICES PROVIDED HEREUNDER.
//
// Distributions of NIST software should also include copyright and licensing
// statements of any third-party software that are legally bundled with the
// code in compliance with the conditions of those licenses.
// 
// ================================================================

// ================================================================
// 
// Authors: Derek Juba <derek.juba@nist.gov>
// Created: Fri Feb 13 13:31:22 2015 EDT
//
// ================================================================

#ifndef WALKER_EXTERIOR_H
#define WALKER_EXTERIOR_H

#include <vector>
#include <cmath>

#include "../Geometry/Sphere.h"
#include "../Geometry/Vector3.h"

namespace zeno {

/// Performs random walks starting on a bounding sphere and determines whether 
/// they hit an object, allowing for a given relative error in distance.
///
template <class T, 
  class RandomNumberGenerator,
  class NearestSurfacePointFinder,
  class RandomSpherePointGenerator,
  class BiasedSpherePointGenerator>
class WalkerExterior {
 public:
  WalkerExterior(RandomNumberGenerator * randomNumberGenerator, 
		 Sphere<T> const & boundingSphere, 
		 NearestSurfacePointFinder const & nearestSurfacePointFinder,
		 T shellThickness);

  ~WalkerExterior();

  void walk(bool * hitObject, int * numSteps,
	    Vector3<T> * startPoint, Vector3<T> * endPoint);

 private:
  RandomNumberGenerator * randomNumberGenerator;
  Sphere<T> const * boundingSphere;
  NearestSurfacePointFinder const * nearestSurfacePointFinder;
  T shellThickness;
};

template <class T, 
  class RandomNumberGenerator,
  class NearestSurfacePointFinder,
  class RandomSpherePointGenerator,
  class BiasedSpherePointGenerator>
WalkerExterior<T, 
               RandomNumberGenerator,
               NearestSurfacePointFinder,
               RandomSpherePointGenerator,
               BiasedSpherePointGenerator>::
  WalkerExterior(RandomNumberGenerator * randomNumberGenerator, 
		 Sphere<T> const & boundingSphere, 
		 NearestSurfacePointFinder const & nearestSurfacePointFinder,
		 T shellThickness) :
  randomNumberGenerator(randomNumberGenerator), 
  boundingSphere(&boundingSphere),
  nearestSurfacePointFinder(&nearestSurfacePointFinder),
  shellThickness(shellThickness) {

}

template <class T, 
  class RandomNumberGenerator,
  class NearestSurfacePointFinder,
  class RandomSpherePointGenerator,
  class BiasedSpherePointGenerator>
WalkerExterior<T, 
               RandomNumberGenerator,
               NearestSurfacePointFinder,
               RandomSpherePointGenerator,
               BiasedSpherePointGenerator>::
  ~WalkerExterior() {

}

/// Perform a random walk and determine whether it hits the object, the number
/// of steps it took, its start and end points, and the surface normal of its
/// hit point.
///
template <class T, 
  class RandomNumberGenerator,
  class NearestSurfacePointFinder,
  class RandomSpherePointGenerator,
  class BiasedSpherePointGenerator>
void 
WalkerExterior<T, 
               RandomNumberGenerator,
               NearestSurfacePointFinder,
               RandomSpherePointGenerator,
               BiasedSpherePointGenerator>::
  walk(bool * hitObject, int * numSteps,
       Vector3<T> * startPoint, Vector3<T> * endPoint) {

  *hitObject = false;
  *numSteps  = 0;

  Vector3<T> position = 
    RandomSpherePointGenerator::generate(randomNumberGenerator, 
					 *boundingSphere);

  *startPoint = position;

  for (;;) {

    T minDistanceSqr = 0;

    nearestSurfacePointFinder->findNearestPoint(position,
						&minDistanceSqr);

    T minDistance = std::sqrt(minDistanceSqr);

    if (minDistance < shellThickness) {
      //walker is absorbed
      
      *endPoint  = position;
      *hitObject = true;
      
      return;
    }
	  
    minDistance += shellThickness;

    (*numSteps)++;

    Sphere<T> stepSphere(position, minDistance);

    position = RandomSpherePointGenerator::generate(randomNumberGenerator, 
						    stepSphere);

    T centerDistSqr = 
      (position - boundingSphere->getCenter()).getMagnitudeSqr();

    if (centerDistSqr > boundingSphere->getRadiusSqr()) {
      //walker left bounding sphere

      T alpha = boundingSphere->getRadius() / std::sqrt(centerDistSqr);

      if (randomNumberGenerator->getRandInRange(0, 1) > (1 - alpha)) {
	//walker is replaced

	position = BiasedSpherePointGenerator::generate(randomNumberGenerator, 
							*boundingSphere,
							position,
							alpha);
      }
      else {
	//walker escapes

	*hitObject = false;
	return;
      }
    }
  }
}

}

#endif

