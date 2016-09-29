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
// Time-stamp: <2016-09-20 17:22:53 dcj>
//
// ================================================================

#ifndef WALKER_EXTERIOR_H
#define WALKER_EXTERIOR_H

#include <vector>

#include "../Geometry/Sphere.h"
#include "../Geometry/Vector3.h"

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
		 T fracErrorBound,
		 T shellThickness);

  ~WalkerExterior();

  void walk(bool * hitObject, int * numSteps,
	    Vector3<T> * startPoint, Vector3<T> * endPoint, 
	    Vector3<T> * normal);

 private:
  RandomNumberGenerator * randomNumberGenerator;
  Sphere<T> const * boundingSphere;
  NearestSurfacePointFinder const * nearestSurfacePointFinder;
  T fracErrorBound;
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
		 T fracErrorBound,
		 T shellThickness) :
  randomNumberGenerator(randomNumberGenerator), 
  boundingSphere(&boundingSphere),
  nearestSurfacePointFinder(&nearestSurfacePointFinder), 
  fracErrorBound(fracErrorBound),
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
       Vector3<T> * startPoint, Vector3<T> * endPoint,
       Vector3<T> * normal) {

  *hitObject = false;
  *numSteps  = 0;

  Vector3<T> position = 
    RandomSpherePointGenerator::generate(randomNumberGenerator, 
					 *boundingSphere);

  *startPoint = position;

  for (;;) {

    T minDistance = 0;
  
    nearestSurfacePointFinder->findNearestPoint(position,
						fracErrorBound, 
						normal,
						&minDistance);

    if (minDistance < shellThickness) {
      //walker is absorbed

      *endPoint  = position;
      *hitObject = true;
      return;
    }

    (*numSteps)++;

    Sphere<T> stepSphere(position, minDistance);

    position = RandomSpherePointGenerator::generate(randomNumberGenerator, 
						    stepSphere);

    T centerDistSqr = 
      (position - boundingSphere->getCenter()).getMagnitudeSqr();

    if (centerDistSqr > boundingSphere->getRadiusSqr()) {
      //walker left bounding sphere

      T alpha = boundingSphere->getRadius() / sqrt(centerDistSqr);

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

#endif

// ================================================================

// Local Variables:
// time-stamp-line-limit: 30
// End:
