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
// Created: Thu Feb 12 17:24:57 2015 EDT
//
// ================================================================

#ifndef BIASED_SPHERE_POINT_DIRECT_H
#define BIASED_SPHERE_POINT_DIRECT_H

#include <cmath>
#include <cassert>

#include "../Geometry/Sphere.h"
#include "../Geometry/Vector3.h"

namespace zeno {

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

  T sinTheta = std::sqrt(1 - std::pow(cosTheta, 2));

  T phi = rng->getRandInRange(0, 2*M_PI);

  Vector3<T> point(sinTheta * std::cos(phi),
		   sinTheta * std::sin(phi),
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

  
  T num = (-std::pow(1 - alpha, 2) + 
	   2*(1 - alpha)*(1 + std::pow(alpha, 2))*R + 
	   2*alpha*(1 + std::pow(alpha, 2))*std::pow(R, 2));

  T den = std::pow(1 - alpha + 2*alpha*R, 2);

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

}

#endif

