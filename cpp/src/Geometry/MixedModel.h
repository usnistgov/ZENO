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
// Created: 2017-04-05
//
// ================================================================

#ifndef MIXEDMODEL_H
#define MIXEDMODEL_H

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "Sphere.h"
#include "Cuboid.h"
#include "Voxels.h"
#include "Vector3.h"
#include "AABBTree.h"

#include <vector>
#include <string>
#include <limits>

// ================================================================

/// Represents a geometric model consisting of a mixture of different types of
/// geometric primitives.
///
template <class T>
class MixedModel {
  using SphereT   = Sphere<T>;
  using CuboidT   = Cuboid<T>;
  using PointT    = Vector3<T>;
  
 public:
  MixedModel();
  ~MixedModel();

  void addSphere(SphereT const & sphere);
  void addCuboid(CuboidT const & cuboid);

  bool loadVoxels(std::string const & inputFileName);

  void mpiSend() const;
  void mpiReceive();

  bool isEmpty() const;
  
  void preprocess();

  void findNearestPoint(PointT const & queryPoint,
			T * minDistanceSqr) const;

  bool contains(PointT const & queryPoint) const;

  PointT findFarthestPoint(PointT const & queryPoint) const;

  CuboidT generateAABB() const;
  
  std::vector<SphereT> const * getSpheres() const;
  std::vector<CuboidT> const * getCuboids() const;
  
 private:
  template <class Primitive>
  void getFarthestPoint(std::vector<Primitive> const & primitives,
			PointT const & queryPoint,
			PointT * farthestPoint) const;

  template <class Primitive>
  void generateAABB(std::vector<Primitive> const & primitives,
                    CuboidT * aabb) const;
    
  AABBTree<SphereT> aabbTreeSpheres;
  AABBTree<CuboidT> aabbTreeCuboids;

  std::vector<SphereT> spheres;
  std::vector<CuboidT> cuboids;

  bool preprocessed;

  void mpiBcastLargeArray(long long int arraySize,
			  double * array) const;
};

// ================================================================

/// Construct an empty model.
///
template <class T>
MixedModel<T>::MixedModel()
: aabbTreeSpheres(),
  aabbTreeCuboids(),
  spheres(),
  cuboids(),
  preprocessed(false) {

}

template <class T>
MixedModel<T>::~MixedModel() {

}

/// Add a sphere to the model.
///
template <class T>
void
MixedModel<T>::addSphere(SphereT const & sphere) {
  if (sphere.getRadius() == 0) {
    std::cerr << "Warning: degenerate sphere" << std::endl
	      << sphere << std::endl;
  }
  
  spheres.push_back(sphere);

  preprocessed = false;
}

/// Add a cuboid to the model.
///
template <class T>
void
MixedModel<T>::addCuboid(CuboidT const & cuboid) {
  if (cuboid.getMaxCoords() == cuboid.getMinCoords()) {
    std::cerr << "Warning: degenerate cuboid" << std::endl
	      << cuboid << std::endl;
  }
  
  cuboids.push_back(cuboid);

  preprocessed = false;
}

/// Add cuboids to the model representing voxels stored in a .fits.gz file.
///
template <class T>
bool
MixedModel<T>::loadVoxels(std::string const & inputFileName) {
  Voxels<T> voxels;

  bool loadSuccess = voxels.loadFitsGz(inputFileName);

  if (!loadSuccess) {
    return false;
  }
  
  voxels.generateCuboids(&cuboids);

  preprocessed = false;
  
  return loadSuccess;
}

/// Broadcasts the geometric primitives stored in the model over MPI.
///
template <class T>
void
MixedModel<T>::mpiSend() const {
#ifdef USE_MPI
  long long int sizeArray[2];

  sizeArray[0] = spheres.size();
  sizeArray[1] = cuboids.size();
  
  MPI_Bcast(sizeArray, 2, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

  long long int geometryArraySize =
    spheres.size() * 4 +
    cuboids.size() * 6;
  
  double * geometryArray = new double[geometryArraySize];

  long long int geometryArrayIndex = 0;

  // spheres
  for (typename std::vector<SphereT>::const_iterator it = spheres.begin();
       it != spheres.end();
       ++it) {

    for (int dim = 0; dim < 3; ++ dim) {
      geometryArray[geometryArrayIndex++] = it->getCenter().get(dim);
    }
    
    geometryArray[geometryArrayIndex++] = it->getRadius();
  }

  // cuboids
  for (typename std::vector<CuboidT>::const_iterator it = cuboids.begin();
       it != cuboids.end();
       ++it) {

    for (int dim = 0; dim < 3; ++ dim) {
      geometryArray[geometryArrayIndex++] = it->getMinCoords().get(dim);
    }

    for (int dim = 0; dim < 3; ++ dim) {
      geometryArray[geometryArrayIndex++] = it->getMaxCoords().get(dim);
    }
  }
  
  assert(geometryArrayIndex == geometryArraySize);

  mpiBcastLargeArray(geometryArraySize, geometryArray);

  delete [] geometryArray;
#endif
}

/// Receives a set of geometric primitives over MPI and stores them in the
/// model.
///
template <class T>
void
MixedModel<T>::mpiReceive() {
#ifdef USE_MPI
  long long int sizeArray[2];

  MPI_Bcast(&sizeArray, 2, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

  long long int spheresSize   = sizeArray[0];
  long long int cuboidsSize   = sizeArray[1];
  
  long long int geometryArraySize =
    spheresSize * 4 +
    cuboidsSize * 6;

  double * geometryArray = new double[geometryArraySize];

  mpiBcastLargeArray(geometryArraySize, geometryArray);

  spheres.clear();
  cuboids.clear();

  spheres.reserve(spheresSize);
  cuboids.reserve(cuboidsSize);
  
  long long int geometryArrayIndex = 0;

  // spheres
  for (long long int i = 0; i < spheresSize; i++) {
    Vector3<T> center;
    
    for (int dim = 0; dim < 3; ++ dim) {
      center.set(dim, geometryArray[geometryArrayIndex++]);
    }

    T radius = geometryArray[geometryArrayIndex++];
    
    spheres.emplace_back(center, radius);
  }

  // cuboids
  for (long long int i = 0; i < cuboidsSize; i++) {
    Vector3<T> minCoords;

    for (int dim = 0; dim < 3; ++ dim) {
      minCoords.set(dim, geometryArray[geometryArrayIndex++]);
    }

    Vector3<T> maxCoords;

    for (int dim = 0; dim < 3; ++ dim) {
      maxCoords.set(dim, geometryArray[geometryArrayIndex++]);
    }

    cuboids.emplace_back(minCoords, maxCoords);
  }

  assert(geometryArrayIndex == geometryArraySize);

  delete [] geometryArray;

  preprocessed = false;
#endif
}

/// Returns whether the model contains no geometric primitives.
///
template <class T>
bool
MixedModel<T>::isEmpty() const {
  return (spheres.empty() &&
	  cuboids.empty());
}

/// Preprocesses geometric primitives stored in the model for use with
/// searching functions.
///
template <class T>
void
MixedModel<T>::preprocess() {
  aabbTreeSpheres.preprocess(&spheres);
  aabbTreeCuboids.preprocess(&cuboids);

  preprocessed = true;
}

/// Finds the nearest point on the model to the query point.  Computes the
/// distance squared to the point.
///
template <class T>
void
MixedModel<T>::findNearestPoint(PointT const & queryPoint,
				T * minDistanceSqr) const {
  assert(preprocessed);

  *minDistanceSqr = std::numeric_limits<T>::infinity();
  
  SphereT const * nearestSphere = NULL;
  
  aabbTreeSpheres.findNearestObject(queryPoint,
 				    &nearestSphere, 
				    minDistanceSqr);

  CuboidT const * nearestCuboid = NULL;

  aabbTreeCuboids.findNearestObject(queryPoint,
				    &nearestCuboid, 
				    minDistanceSqr);
}

/// Returns whether the model contains the query point.
///
template <class T>
bool
MixedModel<T>::contains(PointT const & queryPoint) const {
  assert(preprocessed);

  return (aabbTreeSpheres.objectsContain(queryPoint) ||
	  aabbTreeCuboids.objectsContain(queryPoint));
}

/// Returns the point on the surface of the model that is farthest from the
/// query point.
///
template <class T>
typename MixedModel<T>::PointT
MixedModel<T>::findFarthestPoint(PointT const & queryPoint) const {

  PointT farthestPoint = queryPoint;

  getFarthestPoint(spheres,
		   queryPoint,
		   &farthestPoint);

  getFarthestPoint(cuboids,
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
MixedModel<T>::getFarthestPoint(std::vector<Primitive> const & primitives,
				PointT const & queryPoint,
				PointT * farthestPoint) const {

  for (typename std::vector<Primitive>::const_iterator it = primitives.begin();
       it != primitives.end();
       ++it) {

    PointT primitiveFarthestPoint = it->getFarthestPoint(queryPoint);

    if ((primitiveFarthestPoint - queryPoint).getMagnitudeSqr() >
	(*farthestPoint - queryPoint).getMagnitudeSqr()) {

      *farthestPoint = primitiveFarthestPoint;
    }
  }
}

/// Returns the tightest Axis Aligned Bounding Box of the model.
///
template <class T>
typename MixedModel<T>::CuboidT
MixedModel<T>::generateAABB() const {

  CuboidT aabb(Vector3<T>(std::numeric_limits<T>::max(), 
		 	  std::numeric_limits<T>::max(),
			  std::numeric_limits<T>::max()),

	       Vector3<T>(std::numeric_limits<T>::lowest(),
			  std::numeric_limits<T>::lowest(),
			  std::numeric_limits<T>::lowest()));

  generateAABB(spheres,
	       &aabb);

  generateAABB(cuboids,
	       &aabb);

  return aabb;
}

/// Sets aabb to the tightest Axis Aligned Bounding Box of the given
/// primitives and the given initial aabb.
///
template <class T>
template <class Primitive>
void 
MixedModel<T>::generateAABB(std::vector<Primitive> const & primitives,
                            CuboidT * aabb) const {

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
std::vector<typename MixedModel<T>::SphereT> const *
MixedModel<T>::getSpheres() const {
  return &spheres;
}

template <class T>
std::vector<typename MixedModel<T>::CuboidT> const *
MixedModel<T>::getCuboids() const {
  return &cuboids;
}

/// Uses MPI_Bcast to send or receive an array with size represented by
/// long long int instead of int
///
template <class T>
void
MixedModel<T>::mpiBcastLargeArray(long long int arraySize,
				  double * array) const {
#ifdef USE_MPI
  const long long int maxBroadcastCount = std::numeric_limits<int>::max();

  long long int countBroadcastSoFar = 0;

  while (countBroadcastSoFar < arraySize) {
    long long int countLeftToBroadcast =
      arraySize - countBroadcastSoFar;

    int countToBroadcast =
      (int)std::min(countLeftToBroadcast, maxBroadcastCount);

    double * arrayLeftToBroadcast = array + countBroadcastSoFar;
    
    MPI_Bcast(arrayLeftToBroadcast, countToBroadcast,
	      MPI_DOUBLE, 0, MPI_COMM_WORLD);

    countBroadcastSoFar += countToBroadcast;
  }
#endif
}

#endif
