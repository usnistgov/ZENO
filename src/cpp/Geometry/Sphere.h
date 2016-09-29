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
// Date:    Thu Nov 13 12:27:37 2014 EDT
//
// Time-stamp: <2016-09-19 15:30:56 dcj>
//
// ================================================================

#ifndef SPHERE_H
#define SPHERE_H

#include <cmath>
#include <ostream>

#include "Vector3.h"

// ================================================================

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
Sphere<T>::Sphere() {

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
  return center.get(dim) + radius;
}

/// Returns the minimum coordinate value of any point on the sphere along the
/// given dimension.
///
template <class T>
T 
Sphere<T>::getMinCoord(int dim) const {
  return center.get(dim) - radius;
}

template <class T>
T 
Sphere<T>::getVolume() const {
  return (4./3.)*M_PI*pow(radius, 3);
}

// ================================================================

#endif

// ================================================================

// Local Variables:
// time-stamp-line-limit: 30
// End:
