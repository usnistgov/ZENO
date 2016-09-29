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
// Date:    Thu Feb 12 17:24:57 2015 EDT
//
// Time-stamp: <2016-09-19 17:43:13 dcj>
//
// ================================================================

#ifndef RANDOM_SPHERE_POINT_MARSAGLIA_H
#define RANDOM_SPHERE_POINT_MARSAGLIA_H

#include "../Geometry/Sphere.h"
#include "../Geometry/Vector3.h"

/// Generates random sample points on a sphere from a uniform
/// distribution using the Marsaglia method. 
///
/// Marsaglia, George. Choosing a Point from the Surface of a Sphere. Ann. Math. Statist. 43 (1972), no. 2, 645--646. doi:10.1214/aoms/1177692644. http://projecteuclid.org/euclid.aoms/1177692644.
///
template <class T, class RNG>
class RandomSpherePointMarsaglia {
 public:
  static Vector3<T> generate(RNG * rng, Sphere<T> const & sphere);
};

/// Generates a random point on the given sphere.
///
template <class T, class RNG>
Vector3<T> 
RandomSpherePointMarsaglia<T, RNG>::generate(RNG * rng, 
					     Sphere<T> const & sphere) {

  //generate point uniformly distributed on unit sphere

  T x = 0, y = 0, lengthSqr = 0;

  do {
    x = rng->getRandInRange(-1, 1);
    y = rng->getRandInRange(-1, 1);
    
    lengthSqr = x*x + y*y;
  } 
  while (lengthSqr > 1);

  T scale = 2*sqrt(1 - lengthSqr);

  Vector3<T> spherePoint(scale*x,
			 scale*y,
			 2*lengthSqr - 1);

  //move point onto requested sphere

  spherePoint *= sphere.getRadius();

  spherePoint += sphere.getCenter();

  return spherePoint;
}

#endif

// ================================================================

// Local Variables:
// time-stamp-line-limit: 30
// End:
