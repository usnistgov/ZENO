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
// Created: Thu Nov 13 12:27:37 2014 EDT
//
// ================================================================

#ifndef CUBOID_H
#define CUBOID_H

#include <cmath>
#include <ostream>
#include <algorithm>

#include "Vector3.h"

// ================================================================

/// Represents a cuboid.
///
template <class T>
class Cuboid {
public:
  Cuboid();
  Cuboid(Cuboid<T> const & b);
  Cuboid(Vector3<T> const & minCoords, Vector3<T> const & maxCoords);

  void setMaxCoords(Vector3<T> const & maxCoords);
  Vector3<T> getMaxCoords() const;
  
  void setMinCoords(Vector3<T> const & minCoords);
  Vector3<T> getMinCoords() const;
 
  void setMaxCoord(int dim, T maxCoord);
  T getMaxCoord(int dim) const;
  
  void setMinCoord(int dim, T minCoord);
  T getMinCoord(int dim) const;

  Vector3<T> getCenter() const;
  
  T getVolume() const;

  T getDistanceTo(Vector3<T> const & point) const;
  T getDistanceSqrTo(Vector3<T> const & point) const;

  Vector3<T> getFarthestPoint(Vector3<T> const & queryPoint) const;
  
  bool contains(Vector3<T> const & point) const;
  
private:
  Vector3<T> maxCoords;
  Vector3<T> minCoords;
};

template <class T>
std::ostream & operator<<(std::ostream & os, const Cuboid<T> & rhs)
{
  os << rhs.getMinCoords() << ", " << rhs.getMaxCoords();

  return os;
}

template <class T>
Cuboid<T>::Cuboid()
: maxCoords(),
  minCoords() {

}

template <class T>
Cuboid<T>::Cuboid(Cuboid<T> const & b) 
: maxCoords(b.maxCoords),
  minCoords(b.minCoords) {

}

template <class T>
Cuboid<T>::Cuboid(Vector3<T> const & minCoords, Vector3<T> const & maxCoords) 
: maxCoords(maxCoords),
  minCoords(minCoords) {

}

template <class T>
void
Cuboid<T>::setMaxCoords(Vector3<T> const & maxCoords) {
  this->maxCoords = maxCoords;
}

template <class T>
Vector3<T>
Cuboid<T>::getMaxCoords() const {
  return maxCoords;
}

template <class T>
void
Cuboid<T>::setMinCoords(Vector3<T> const & minCoords) {
  this->minCoords = minCoords;
}

template <class T>
Vector3<T>
Cuboid<T>::getMinCoords() const {
  return minCoords;
}

/// Sets the maximum coordinate value of the cuboid along the
/// given dimension.
///
template <class T>
void 
Cuboid<T>::setMaxCoord(int dim, T maxCoord) {
  maxCoords[dim] = maxCoord;
}

/// Returns the maximum coordinate value the cuboid along the
/// given dimension.
///
template <class T>
T 
Cuboid<T>::getMaxCoord(int dim) const {
  return maxCoords[dim];
}

/// Sets the minimum coordinate value of the cuboid along the
/// given dimension.
///
template <class T>
void
Cuboid<T>::setMinCoord(int dim, T minCoord) {
  minCoords[dim] = minCoord;
}

/// Returns the minimum coordinate value of the cuboid along the
/// given dimension.
///
template <class T>
T 
Cuboid<T>::getMinCoord(int dim) const {
  return minCoords[dim];
}

template <class T>
Vector3<T>
Cuboid<T>::getCenter() const {
  return (maxCoords + minCoords) / 2;
}

template <class T>
T 
Cuboid<T>::getVolume() const {
  return
    (maxCoords.getX() - minCoords.getX()) *
    (maxCoords.getY() - minCoords.getY()) *
    (maxCoords.getZ() - minCoords.getZ());
}

/// Returns the distance from the point to the surface of the cuboid.
/// Returns 0 if the point is inside the cuboid.
///
template <class T>
T
Cuboid<T>::getDistanceTo(Vector3<T> const & point) const {
  return sqrt(getDistanceSqrTo(point));
}

/// Returns the distance squared from the point to the surface of the cuboid.
/// Returns 0 if the point is inside the cuboid.
///
template <class T>
T
Cuboid<T>::getDistanceSqrTo(Vector3<T> const & point) const {
  T distanceSqr = 0;

  for (int dim = 0; dim < 3; ++dim) {
    T coordDistance = std::max(minCoords[dim] - point[dim],
			       point[dim] - maxCoords[dim]);

    coordDistance = std::max(coordDistance, (T)0);

    distanceSqr += pow(coordDistance, 2);
  }

  return distanceSqr;
}

/// Returns the point on the surface of the cuboid that is farthest from the
/// given query point.
///
template <class T>
Vector3<T>
Cuboid<T>::getFarthestPoint(Vector3<T> const & queryPoint) const {
  Vector3<T> farthestPoint;

  Vector3<T> center = (minCoords + maxCoords) * 0.5;
  
  for (int dim = 0; dim < 3; ++ dim) {
    if (queryPoint[dim] < center[dim]) {
      farthestPoint[dim] = maxCoords[dim];
    }
    else {
      farthestPoint[dim] = minCoords[dim];
    }
  }

  return farthestPoint;
}

/// Returns whether the cuboid encloses the point.
///
template <class T>
bool
Cuboid<T>::contains(Vector3<T> const & point) const {
  for (int dim = 0; dim < 3; ++dim) {
    if (minCoords[dim] > point[dim] ||
	maxCoords[dim] < point[dim]) {
      
      return false;
    }
  }

  return true;
}

// ================================================================

#endif

