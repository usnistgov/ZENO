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
// Time-stamp: <2016-09-19 15:06:06 dcj>
//
// ================================================================

#ifndef RANDOM_BALL_POINT_REJECTION_H
#define RANDOM_BALL_POINT_REJECTION_H

#include "../Geometry/Sphere.h"
#include "../Geometry/Vector3.h"

/// Generates random sample points inside a sphere from a uniform
/// distribution using a rejection (iterative) method. 
///
template <class T, class RNG>
class RandomBallPointRejection {
 public:
  static Vector3<T> generate(RNG * rng, Sphere<T> const & ball);
};

/// Generates a random point inside the given sphere.
///
template <class T, class RNG>
Vector3<T> 
RandomBallPointRejection<T, RNG>::generate(RNG * rng, 
					   Sphere<T> const & ball) {

  //generate point uniformly distributed on unit ball

  T x = 0, y = 0, z = 0, lengthSqr = 0;

  do {
    x = rng->getRandInRange(-1, 1);
    y = rng->getRandInRange(-1, 1);
    z = rng->getRandInRange(-1, 1);
    
    lengthSqr = x*x + y*y + z*z;
  } 
  while (lengthSqr > 1);

  Vector3<T> ballPoint(x, y, z);

  //move point onto requested ball

  ballPoint *= ball.getRadius();

  ballPoint += ball.getCenter();

  return ballPoint;
}

#endif

// ================================================================

// Local Variables:
// time-stamp-line-limit: 30
// End:
