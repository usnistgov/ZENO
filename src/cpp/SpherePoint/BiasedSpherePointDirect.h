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
// Time-stamp: <2016-09-19 14:57:03 dcj>
//
// ================================================================

#ifndef BIASED_SPHERE_POINT_DIRECT_H
#define BIASED_SPHERE_POINT_DIRECT_H

#include <cmath>
#include <cassert>

#include "../Geometry/Sphere.h"
#include "../Geometry/Vector3.h"

/// Generates random sample points on a sphere from a biased (non-uniform)
/// distribution using a direct (non-iterative) method.  The distribution is
/// suitable for reinserting random walkers that have left the launch sphere. 
///
template <class T, 
          class RandomNumberGenerator, 
          class RandomSpherePointGenerator>
class BiasedSpherePointDirect {
 public:
  static Vector3<T> generate(RandomNumberGenerator * rng, 
			     Sphere<T> const & sphere,
			     Vector3<T> const & distributionCenter,
			     T alpha);

 private:
  static T computeCosTheta(T alpha, T R);
};

/// Generates a random point on the given sphere from a distribution centered
/// at the given point with distribution parameter "alpha".
///
template <class T, 
          class RandomNumberGenerator, 
          class RandomSpherePointGenerator>
Vector3<T> 
BiasedSpherePointDirect<T, 
                        RandomNumberGenerator,
                        RandomSpherePointGenerator>::
  generate(RandomNumberGenerator * rng, 
	   Sphere<T> const & sphere,
	   Vector3<T> const & distributionCenter,
	   T alpha) {

  //generate point on unit sphere centered at origin distributed around Z axis

  T R = rng->getRandInRange(0, 1);

  T cosTheta = computeCosTheta(alpha, R);

  T sinTheta = sqrt(1 - pow(cosTheta, 2));

  T phi = rng->getRandInRange(0, 2*M_PI);

  Vector3<T> point(sinTheta * cos(phi),
		   sinTheta * sin(phi),
		   cosTheta);

  //rotate point to be distributed around requested distribution center

  Vector3<T> recenteredDistributionCenter =
    distributionCenter - sphere.getCenter();

  recenteredDistributionCenter.normalize();

  //will cause problems if recenteredDistributionCenter is close to <0, 0, 1>
  Vector3<T> rotationAxis =
    recenteredDistributionCenter.cross(Vector3<T>(0, 0, 1));

  T sinRotationAngle = rotationAxis.getMagnitude();

  T cosRotationAngle =
    recenteredDistributionCenter.dot(Vector3<T>(0, 0, 1));

  point.rotate(cosRotationAngle, sinRotationAngle, rotationAxis);

  //move point onto requested sphere

  point *= sphere.getRadius();
  point += sphere.getCenter();

  return point;
}

template <class T, 
          class RandomNumberGenerator, 
          class RandomSpherePointGenerator>
T 
BiasedSpherePointDirect<T, 
                        RandomNumberGenerator,
                        RandomSpherePointGenerator>::
  computeCosTheta(T alpha, T R) {

  
    T num = (-pow(1 - alpha, 2) + 
	     2*(1 - alpha)*(1 + pow(alpha, 2))*R + 
	     2*alpha*(1 + pow(alpha, 2))*pow(R, 2));

    T den = pow(1 - alpha + 2*alpha*R, 2);

    assert(den != 0);

    T cosTheta = num/den;

    //correct for numerical errors
    if (cosTheta > 1) {
      // std::cerr << "Warning: cosTheta > 1 (cosTheta == " << cosTheta << ")" 
      //           << std::endl;

      cosTheta = 1;
    }

    return cosTheta;
}

#endif

// ================================================================

// Local Variables:
// time-stamp-line-limit: 30
// End:
