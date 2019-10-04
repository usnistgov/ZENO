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

#ifndef SPHERE_H
#define SPHERE_H

#include <cmath>
#include <ostream>

#include "Vector3.h"

// ================================================================

namespace zeno {
  
/// Represents a sphere.
///
template <class T>
class Sphere {
public:
  Sphere();
  Sphere(Sphere<T> const & b);
  Sphere(Vector3<T> const & center, T radius);

  void setCenter(Vector3<T> const & center);
  Vector3<T> getCenter() const;

  void setRadiusSqr(T radiusSqr);
  T getRadiusSqr() const;

  void setRadius(T radius);
  T getRadius() const;

  T getMaxCoord(int dim) const;
  T getMinCoord(int dim) const;

  T getVolume() const;

  T getDistanceTo(Vector3<T> const & point) const;
  T getDistanceSqrTo(Vector3<T> const & point) const;

  Vector3<T> getFarthestPoint(Vector3<T> const & queryPoint) const;
  
  bool contains(Vector3<T> const & point) const;
  bool contains(Sphere<T> const & sphere) const;
  
private:
  Vector3<T> center;
  T radius;
};

template <class T>
std::ostream & operator<<(std::ostream & os, const Sphere<T> & rhs)
{
  os << rhs.getCenter() << ", " << rhs.getRadius();

  return os;
}

template <class T>
Sphere<T>::Sphere()
: center(),
  radius() {

}

template <class T>
Sphere<T>::Sphere(Sphere<T> const & b) 
: center(b.center),
  radius(b.radius) {

}

template <class T>
Sphere<T>::Sphere(Vector3<T> const & center, T radius) 
: center(center),
  radius(radius) {

}

template <class T>
void 
Sphere<T>::setCenter(Vector3<T> const & center) {
  this->center = center;
}

template <class T>
Vector3<T> 
Sphere<T>::getCenter() const {
  return center;
}

template <class T>
void 
Sphere<T>::setRadiusSqr(T radiusSqr) {
  this->radius = sqrt(radiusSqr);
}

template <class T>
T 
Sphere<T>::getRadiusSqr() const {
  return pow(radius, 2);
}

template <class T>
void 
Sphere<T>::setRadius(T radius) {
  this->radius = radius;
}

template <class T>
T 
Sphere<T>::getRadius() const {
  return radius;
}

/// Returns the maximum coordinate value of any point on the sphere along the
/// given dimension.
///
template <class T>
T 
Sphere<T>::getMaxCoord(int dim) const {
  return center[dim] + radius;
}

/// Returns the minimum coordinate value of any point on the sphere along the
/// given dimension.
///
template <class T>
T 
Sphere<T>::getMinCoord(int dim) const {
  return center[dim] - radius;
}

template <class T>
T 
Sphere<T>::getVolume() const {
  return (4./3.)*M_PI*pow(radius, 3);
}

/// Returns the distance from the point to the surface of the sphere.
/// Distances are positive for points outside the sphere and negative for
/// points inside the sphere.
///
template <class T>
T
Sphere<T>::getDistanceTo(Vector3<T> const & point) const {
  return (point - center).getMagnitude() - radius;
}

/// Returns the distance squared from the point to the surface of the sphere.
///
template <class T>
T
Sphere<T>::getDistanceSqrTo(Vector3<T> const & point) const {
  return pow(getDistanceTo(point), 2);
}

/// Returns the point on the surface of the sphere that is farthest from the
/// given query point.
///
template <class T>
Vector3<T>
Sphere<T>::getFarthestPoint(Vector3<T> const & queryPoint) const {
  // If center == queryPoint then there are many farthest points, so return
  // one arbitrarily
  if (center == queryPoint) {
    return center + Vector3<T>(radius, 0, 0);
  }
  
  Vector3<T> radiusVec = (center - queryPoint).normalized() * radius;

  Vector3<T> farthestPoint = center + radiusVec; 

  return farthestPoint;
}

/// Returns whether the sphere encloses the point.
///
template <class T>
bool
Sphere<T>::contains(Vector3<T> const & point) const {
  return (point - center).getMagnitudeSqr() < pow(radius, 2);
}

/// Returns whether the sphere completely encloses the given sphere.
///
template <class T>
bool
Sphere<T>::contains(Sphere<T> const & sphere) const {
  return (sphere.center - center).getMagnitude() + sphere.radius <= radius;
}

}

// ================================================================

#endif

