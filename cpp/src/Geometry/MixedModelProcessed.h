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
// Created: 2019-05-20
//
// ================================================================

#ifndef MIXED_MODEL_PROCESSED_H
#define MIXED_MODEL_PROCESSED_H

#include "Sphere.h"
#include "Cuboid.h"
#include "Triangle.h"
#include "Voxels.h"
#include "Vector3.h"
#include "AABBTree.h"
#include "MixedModel.h"

#include <vector>
#include <string>
#include <limits>

// ================================================================

namespace zeno {
  
/// Represents a geometric model consisting of a mixture of different types of
/// geometric primitives.
///
template <class T>
class MixedModelProcessed {
 public:
  MixedModelProcessed();
  ~MixedModelProcessed();

  void addMixedModel(MixedModel<T> * mixedModel);

  bool isEmpty() const;
  
  void preprocess();

  void findNearestPoint(Vector3<T> const & queryPoint,
			T * minDistanceSqr) const;

  bool contains(Vector3<T> const & queryPoint) const;

  Vector3<T> findFarthestPoint(Vector3<T> const & queryPoint) const;

  Cuboid<T> generateAABB() const;
  
  std::vector<Sphere<T> > const * getSpheres() const;
  std::vector<Cuboid<T> > const * getCuboids() const;
  std::vector<Triangle<T> > const * getTriangles() const;
  
 private:
  template <class Primitive>
  void getFarthestPoint(std::vector<Primitive> const & primitives,
			Vector3<T> const & queryPoint,
			Vector3<T> * farthestPoint) const;

  template <class Primitive>
  void generateAABB(std::vector<Primitive> const & primitives,
                    Cuboid<T> * aabb) const;
    
  AABBTree<Sphere<T> > aabbTreeSpheres;
  AABBTree<Cuboid<T> > aabbTreeCuboids;
  AABBTree<Triangle<T> > aabbTreeTriangles;

  std::vector<Sphere<T> > * spheres;
  std::vector<Cuboid<T> > * cuboids;
  std::vector<Triangle<T> > * triangles;

  bool preprocessed;
};

// ================================================================

/// Construct an empty model.
///
template <class T>
MixedModelProcessed<T>::MixedModelProcessed()
: aabbTreeSpheres(),
  aabbTreeCuboids(),
  aabbTreeTriangles(),
  spheres(nullptr),
  cuboids(nullptr),
  triangles(nullptr),
  preprocessed(false) {

}

template <class T>
MixedModelProcessed<T>::~MixedModelProcessed() {

}

template <class T>
void
MixedModelProcessed<T>::addMixedModel(MixedModel<T> * mixedModel) {
  spheres = mixedModel->getAndLockSpheres();
  cuboids = mixedModel->getAndLockCuboids();
  triangles = mixedModel->getAndLockTriangles();

  preprocessed = false;
}

/// Returns whether the model contains no geometric primitives.
///
template <class T>
bool
MixedModelProcessed<T>::isEmpty() const {
  return ((spheres == nullptr || spheres->empty()) &&
	  (cuboids == nullptr || cuboids->empty()) &&
	  (triangles == nullptr || triangles->empty()));
}

/// Preprocesses geometric primitives stored in the model for use with
/// searching functions.
///
template <class T>
void
MixedModelProcessed<T>::preprocess() {
  aabbTreeSpheres.preprocess(spheres);
  aabbTreeCuboids.preprocess(cuboids);
  aabbTreeTriangles.preprocess(triangles);

  preprocessed = true;
}

/// Finds the nearest point on the model to the query point.  Computes the
/// distance squared to the point.
///
template <class T>
void
MixedModelProcessed<T>::findNearestPoint(Vector3<T> const & queryPoint,
					 T * minDistanceSqr) const {
  assert(preprocessed);

  *minDistanceSqr = std::numeric_limits<T>::infinity();
  
  Sphere<T> const * nearestSphere = NULL;
  
  aabbTreeSpheres.findNearestObject(queryPoint,
 				    &nearestSphere, 
				    minDistanceSqr);

  Cuboid<T> const * nearestCuboid = NULL;

  aabbTreeCuboids.findNearestObject(queryPoint,
				    &nearestCuboid, 
				    minDistanceSqr);

  Triangle<T> const * nearestTriangle = NULL;

  aabbTreeTriangles.findNearestObject(queryPoint,
				      &nearestTriangle, 
				      minDistanceSqr);
}

/// Returns whether the model contains the query point.
///
template <class T>
bool
MixedModelProcessed<T>::contains(Vector3<T> const & queryPoint) const {
  assert(preprocessed);

  T minDistanceSqr = std::numeric_limits<T>::infinity();
  
  Triangle<T> const * nearestTriangle = NULL;

  aabbTreeTriangles.findNearestObject(queryPoint,
				      &nearestTriangle, 
				      &minDistanceSqr);

  bool meshContains = ((nearestTriangle != NULL) &&
		       (nearestTriangle->isAbove(queryPoint))); 
  
  return (meshContains ||
	  aabbTreeSpheres.objectsContain(queryPoint) ||
	  aabbTreeCuboids.objectsContain(queryPoint));
}

/// Returns the point on the surface of the model that is farthest from the
/// query point.
///
template <class T>
Vector3<T>
MixedModelProcessed<T>::findFarthestPoint(Vector3<T> const & queryPoint) const {

  Vector3<T> farthestPoint = queryPoint;

  getFarthestPoint(*spheres,
		   queryPoint,
		   &farthestPoint);

  getFarthestPoint(*cuboids,
		   queryPoint,
		   &farthestPoint);

  getFarthestPoint(*triangles,
		   queryPoint,
		   &farthestPoint);

  return farthestPoint;
}

/// Sets farthestPoint to the point on the surface of the given primitives that
/// is farthest from the query point, unless that point is closer to the query
/// point than the initial value of farthestPoint.
///
template <class T>
template <class Primitive>
void
MixedModelProcessed<T>::getFarthestPoint
(std::vector<Primitive> const & primitives,
 Vector3<T> const & queryPoint,
 Vector3<T> * farthestPoint) const {

  for (typename std::vector<Primitive>::const_iterator it = primitives.begin();
       it != primitives.end();
       ++it) {

    Vector3<T> primitiveFarthestPoint = it->getFarthestPoint(queryPoint);

    if ((primitiveFarthestPoint - queryPoint).getMagnitudeSqr() >
	(*farthestPoint - queryPoint).getMagnitudeSqr()) {

      *farthestPoint = primitiveFarthestPoint;
    }
  }
}

/// Returns the tightest Axis Aligned Bounding Box of the model.
///
template <class T>
Cuboid<T>
MixedModelProcessed<T>::generateAABB() const {

  Cuboid<T> aabb(Vector3<T>(std::numeric_limits<T>::max(), 
			    std::numeric_limits<T>::max(),
			    std::numeric_limits<T>::max()),

		 Vector3<T>(std::numeric_limits<T>::lowest(),
			    std::numeric_limits<T>::lowest(),
			    std::numeric_limits<T>::lowest()));

  generateAABB(*spheres,
	       &aabb);

  generateAABB(*cuboids,
	       &aabb);

  generateAABB(*triangles,
	       &aabb);

  return aabb;
}

/// Sets aabb to the tightest Axis Aligned Bounding Box of the given
/// primitives and the given initial aabb.
///
template <class T>
template <class Primitive>
void 
MixedModelProcessed<T>::generateAABB(std::vector<Primitive> const & primitives,
				     Cuboid<T> * aabb) const {

  for (typename std::vector<Primitive>::const_iterator it = primitives.begin();
       it != primitives.end();
       ++it) {

    for (int dim = 0; dim < 3; dim++) {

      T minCoord = it->getMinCoord(dim);
      T maxCoord = it->getMaxCoord(dim);

      if (minCoord < aabb->getMinCoord(dim)) {
	aabb->setMinCoord(dim, minCoord);
      }

      if (maxCoord > aabb->getMaxCoord(dim)) {
	aabb->setMaxCoord(dim, maxCoord);
      }
    }
  }
}

template <class T>
std::vector<Sphere<T> > const *
MixedModelProcessed<T>::getSpheres() const {
  return spheres;
}

template <class T>
std::vector<Cuboid<T> > const *
MixedModelProcessed<T>::getCuboids() const {
  return cuboids;
}

template <class T>
std::vector<Triangle<T> > const *
MixedModelProcessed<T>::getTriangles() const {
  return triangles;
}

}

#endif
