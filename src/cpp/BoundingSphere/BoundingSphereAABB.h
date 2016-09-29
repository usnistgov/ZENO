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
// Date:    Thu Feb 12 16:39:41 2015 EDT
//
// Time-stamp: <2016-09-22 13:15:59 dcj>
//
// ================================================================

#ifndef BOUNDING_SPHERE_AABB_H
#define BOUNDING_SPHERE_AABB_H

#include "../Geometry/Vector3.h"
#include "../Geometry/Sphere.h"

/// Generates the tightest bounding sphere around the tightest
/// Axis-Aligned Bounding-Box around the given object.
///
template <class T>
class BoundingSphereAABB {
 public:
  static Sphere<T> generate(std::vector<Sphere<T> > const & spheres);

 private:
  static Sphere<T> generate(Vector3<T> const & minCornerCoords,
			    Vector3<T> const & maxCornerCoords);
};

/// Returns the bounding sphere for the given set of spheres.
///
template <class T>
Sphere<T> 
BoundingSphereAABB<T>::generate(std::vector<Sphere<T> > const & spheres) {

  Vector3<T> minCornerCoords(std::numeric_limits<T>::max(), 
			     std::numeric_limits<T>::max(),
			     std::numeric_limits<T>::max());

  Vector3<T> maxCornerCoords(0, 0, 0);

  for (typename std::vector<Sphere<T> >::const_iterator it = spheres.begin();
       it != spheres.end();
       ++it) {

    for (int dim = 0; dim < 3; dim++) {

      T minCoord = it->getMinCoord(dim);
      T maxCoord = it->getMaxCoord(dim);

      if (minCoord < minCornerCoords.get(dim)) {
	minCornerCoords.set(dim, minCoord);
      }

      if (maxCoord > maxCornerCoords.get(dim)) {
	maxCornerCoords.set(dim, maxCoord);
      }
    }
  }

  return generate(minCornerCoords,
		  maxCornerCoords);
}

/// Returns the bounding sphere for a box with the given max and min
/// coordinate values.
///
template <class T>
Sphere<T> 
BoundingSphereAABB<T>::generate(Vector3<T> const & minCornerCoords,
				Vector3<T> const & maxCornerCoords) {

  Vector3<T> AABBDiagonal((maxCornerCoords - minCornerCoords) / 2);

  Vector3<T> boundingSphereCenter(minCornerCoords + AABBDiagonal);

  T boundingSphereRadius = AABBDiagonal.getMagnitude();

  Sphere<T> boundingSphere(boundingSphereCenter, boundingSphereRadius);

  return boundingSphere;
}

#endif

// ================================================================

// Local Variables:
// time-stamp-line-limit: 30
// End:
