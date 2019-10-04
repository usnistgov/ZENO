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
// Created: 2017-04-06
//
// ================================================================

#ifndef TRIANGLE_H
#define TRIANGLE_H

#include <cmath>
#include <ostream>
#include <algorithm>

#include "Vector3.h"

// ================================================================

namespace zeno {
  
/// Represents a triangle.
///
template <class T>
class Triangle {
public:
  Triangle();
  
  Triangle(Triangle<T> const & b);
  
  Triangle(Vector3<T> const & v1,
	   Vector3<T> const & v2,
	   Vector3<T> const & v3);

  void setV(Vector3<T> const & v1,
	    Vector3<T> const & v2,
	    Vector3<T> const & v3);
  
  void getV(Vector3<T> * v1,
	    Vector3<T> * v2,
	    Vector3<T> * v3) const;

  Vector3<T> getV1() const;
  Vector3<T> getV2() const;
  Vector3<T> getV3() const;
  
  T getMaxCoord(int dim) const;
  T getMinCoord(int dim) const;

  T getVolume() const;

  T getDistanceTo(Vector3<T> const & point) const;
  T getDistanceSqrTo(Vector3<T> const & point) const;

  Vector3<T> getFarthestPoint(Vector3<T> const & queryPoint) const;
  
  bool contains(Vector3<T> const & point) const;

  bool isAbove(Vector3<T> const & point) const;
  
private:
  Vector3<T> v1, v2, v3;

  Vector3<T> normal;
  
  void preprocess();
};

template <class T>
std::ostream & operator<<(std::ostream & os, const Triangle<T> & rhs)
{
  os << rhs.getV1() << ", " << rhs.getV2() << ", " << rhs.getV3();

  return os;
}

template <class T>
Triangle<T>::Triangle()
: v1(),
  v2(),
  v3(),
  normal() {

}

template <class T>
Triangle<T>::Triangle(Triangle<T> const & b) 
: v1(b.v1),
  v2(b.v2),
  v3(b.v3),
  normal(b.normal) {

}

template <class T>
Triangle<T>::Triangle(Vector3<T> const & v1,
		      Vector3<T> const & v2,
		      Vector3<T> const & v3) 
: v1(v1),
  v2(v2),
  v3(v3),
  normal() {
  
  preprocess();
}

template <class T>
void 
Triangle<T>::setV(Vector3<T> const & v1,
		  Vector3<T> const & v2,
		  Vector3<T> const & v3) {
  this->v1 = v1;
  this->v2 = v2;
  this->v3 = v3;

  preprocess();
}

template <class T>
void 
Triangle<T>::getV(Vector3<T> * v1,
		  Vector3<T> * v2,
		  Vector3<T> * v3) const {
  *v1 = this->v1;
  *v2 = this->v2;
  *v3 = this->v3;
}

template <class T>
Vector3<T> 
Triangle<T>::getV1() const {
  return v1;
}

template <class T>
Vector3<T> 
Triangle<T>::getV2() const {
  return v2;
}

template <class T>
Vector3<T> 
Triangle<T>::getV3() const {
  return v3;
}

template <class T>
void 
Triangle<T>::preprocess() {
  normal = (v2 - v1).cross(v3 - v1);
  normal.normalize();
}

/// Returns the maximum coordinate value of any point on the triangle along the
/// given dimension.
///
template <class T>
T 
Triangle<T>::getMaxCoord(int dim) const {
  return std::max(std::max(v1[dim],
			   v2[dim]),
		  v3[dim]);
}

/// Returns the minimum coordinate value of any point on the triangle along the
/// given dimension.
///
template <class T>
T 
Triangle<T>::getMinCoord(int dim) const {
  return std::min(std::min(v1[dim],
			   v2[dim]),
		  v3[dim]);
}

/// Volume is 0 since a triangle is 2D.
///
template <class T>
T 
Triangle<T>::getVolume() const {
  assert(0);
  
  return 0;
}

/// Returns the distance from the point to the triangle.
///
template <class T>
T
Triangle<T>::getDistanceTo(Vector3<T> const & point) const {
  return sqrt(getDistanceSqrTo(point));
}

/// Returns the distance squared from the point to the triangle.
///
template <class T>
T
Triangle<T>::getDistanceSqrTo(Vector3<T> const & point) const {
  // Eberly, David. "Distance between point and triangle in 3D."
  // Magic Software, (1999).

  Vector3<T> E0 = v2 - v1;
  Vector3<T> E1 = v3 - v1;
  
  Vector3<T> D = v1 - point;

  T a = E0.dot(E0);
  T b = E0.dot(E1);
  T c = E1.dot(E1);
  T d = E0.dot(D);
  T e = E1.dot(D);
  T f = D.dot(D);

  T det = a*c - b*b;
  T s   = b*e - c*d;
  T t   = b*d - a*e;

  if (s + t <= det) {
    if (s < 0) {
      if (t < 0) {
  	// region 4

  	if (d < 0) {
	  t = 0;
	  s = (-d >= a ? 1 : -d/a);
  	}
  	else {
	  s = 0;
	  t = (e >= 0 ? 0 : (-e >= c ? 1 : -e/c));
  	}
      }
      else {
  	// region 3
	
        s = 0;
	t = (e >= 0 ? 0 : (-e >= c ? 1 : -e/c));
      }
    }
    else if (t < 0) {
      // region 5

      t = 0;
      s = (d >= 0 ? 0 : (-d >= a ? 1 : -d/a));
    }
    else {
      // region 0

      T invDet = 1. / det;

      s *= invDet;
      t *= invDet;
    }
  }
  else {
    if (s < 0) {
      // region 2

      T tmp0 = b + d;
      T tmp1 = c + e;

      if (tmp1 > tmp0) {
  	T numer = tmp1 - tmp0;
  	T denom = a - 2*b + c;

	s = (numer >= denom ? 1 : numer/denom);
	t = 1 - s;
      }
      else {
	s = 0;
	t = (tmp1 <= 0 ? 1 : (e >= 0 ? 0 : -e/c));
      }
    }
    else if (t < 0) {
      // region 6
      
      T tmp0 = b + e;
      T tmp1 = a + d;

      if (tmp1 > tmp0) {
  	T numer = tmp1 - tmp0;
  	T denom = a - 2*b + c;

	t = (numer >= denom ? 1 : numer/denom);
	s = 1 - t;
      }
      else {
	t = 0;
	s = (tmp1 <= 0 ? 1 : (d >= 0 ? 0 : -d/a));
      }
    }
    else {
      // region 1

      T numer = c + e - b - d;

      if (numer <= 0) {
  	s = 0;
      }
      else {
  	T denom = a - 2*b + c;

	s = (numer >= denom ? 1 : numer/denom);
      }

      t = 1 - s;
    }
  }

  T distanceSqr = a*s*s + 2*b*s*t + c*t*t + 2*d*s + 2*e*t + f;

  // correct for numerical errors
  if (distanceSqr < 0) {
    //std::cerr << "Warning: distanceSqr < 0 (distanceSqr == " << distanceSqr
    //          << ")" << std::endl;
    
    distanceSqr = 0;
  }

  return distanceSqr;
}

/// Returns the point on the surface of the triangle that is farthest from the
/// given query point.
///
template <class T>
Vector3<T>
Triangle<T>::getFarthestPoint(Vector3<T> const & queryPoint) const {
  T distSqrV1 = (v1 - queryPoint).getMagnitudeSqr();
  T distSqrV2 = (v2 - queryPoint).getMagnitudeSqr();
  T distSqrV3 = (v3 - queryPoint).getMagnitudeSqr();

  Vector3<T> farthestPoint;
  
  if      ((distSqrV1 >= distSqrV2) && (distSqrV1 >= distSqrV3)) {
    farthestPoint = v1;
  }
  else if ((distSqrV2 >= distSqrV1) && (distSqrV2 >= distSqrV3)) {
    farthestPoint = v2;
  }
  else if ((distSqrV3 >= distSqrV1) && (distSqrV3 >= distSqrV2)) {
    farthestPoint = v3;
  }
  else {
    assert(0);
  }

  return farthestPoint;
}

/// Assuming General Position, a triangle will never contain a 3D point
/// since it is 2D.
///
template <class T>
bool
Triangle<T>::contains(Vector3<T> const & point) const {
  assert(0);
  
  return false;
}

/// Returns whether the triangle is above the point.
/// "Above" means that the vector from the point to the trangle's plane is in
/// the same direction as the triangle's normal (which is based on the order of
/// its vertexes).
///
template <class T>
bool
Triangle<T>::isAbove(Vector3<T> const & point) const {
  return ((point - v1).dot(normal) < 0);
}

}

// ================================================================

#endif

