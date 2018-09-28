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
// Created: Thu Feb 12 16:39:41 2015 EDT
//
// ================================================================

#ifndef BOUNDING_SPHERE_GENERATOR_H
#define BOUNDING_SPHERE_GENERATOR_H

#include "Geometry/Vector3.h"
#include "Geometry/Sphere.h"
#include "Geometry/Cuboid.h"
#include "Geometry/MixedModel.h"

/// Generates a tight bounding sphere around the given object.
///
template <class T>
class BoundingSphereGenerator {
 public:
  
  static Sphere<T> generate(MixedModel<T> const & model);

 private:
  
  static Vector3<T> getCenterFromAABB(Cuboid<T> const & aabb);

  static Vector3<T> getCenterRitter(MixedModel<T> const & model);
  
  static T getRadiusFromAABB(Cuboid<T> const & aabb);

  static T getRadiusFromCenter(MixedModel<T> const & model,
			       Vector3<T> const & center);
};

/// Returns the bounding sphere for the given spheres, cuboids, and triangles.
///
template <class T>
Sphere<T> 
BoundingSphereGenerator<T>::generate(MixedModel<T> const & model) {

  Vector3<T> center = getCenterRitter(model);

  T radius = getRadiusFromCenter(model,
		                 center);
  
  return Sphere<T>(center, radius);
}

template <class T>
Vector3<T>
BoundingSphereGenerator<T>::getCenterFromAABB(Cuboid<T> const & aabb) {

  Vector3<T> aabbDiagonal = (aabb.getMaxCoords() -
			     aabb.getMinCoords()) * 0.5;

  Vector3<T> center = aabb.getMinCoords() + aabbDiagonal;

  return center;
}

/// Generates a tight bounding sphere around the given object using Ritter's
/// algorithm.
///
/// Jack Ritter. 1990. An efficient bounding sphere.
/// In Graphics gems, Andrew S. Glassner (Ed.).
/// Academic Press Professional, Inc., San Diego, CA, USA 301-303. 
///
template <class T>
Vector3<T>
BoundingSphereGenerator<T>::getCenterRitter(MixedModel<T> const & model) {

  Cuboid<T> aabb = model.generateAABB();
  
  Vector3<T> initialPoint = getCenterFromAABB(aabb);

  Vector3<T> edgePoint1 = model.findFarthestPoint(initialPoint);

  Vector3<T> edgePoint2 = model.findFarthestPoint(edgePoint1);

  Vector3<T> center = (edgePoint1 + edgePoint2) * 0.5;

  return center;
}

template <class T>
T
BoundingSphereGenerator<T>::getRadiusFromAABB(Cuboid<T> const & aabb) {

  Vector3<T> aabbDiagonal = (aabb.getMaxCoords() -
			     aabb.getMinCoords()) * 0.5;

  T radius = aabbDiagonal.getMagnitude();

  return radius;
}

template <class T>
T
BoundingSphereGenerator<T>::getRadiusFromCenter(MixedModel<T> const & model,
			                        Vector3<T> const & center) {

  Vector3<T> edgePoint = model.findFarthestPoint(center);

  T radius = (edgePoint - center).getMagnitude();

  return radius;
}

#endif

