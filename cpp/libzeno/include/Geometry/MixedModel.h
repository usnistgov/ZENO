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
#include "Triangle.h"
#include "Voxels.h"
#include "Vector3.h"
#include "AABBTree.h"

#include <vector>
#include <string>
#include <limits>

// ================================================================

namespace zeno {
  
/// Represents a geometric model consisting of a mixture of different types of
/// geometric primitives.
///
template <class T>
class MixedModel {
 public:
  MixedModel();
  ~MixedModel();

  void addSphere(Sphere<T> const & sphere);
  void addCuboid(Cuboid<T> const & cuboid);

  bool loadVoxels(std::string const & inputFileName);

  void mpiBroadcast(int root);

  bool isEmpty() const;
  
  std::vector<Sphere<T> > const * getSpheres() const;
  std::vector<Cuboid<T> > const * getCuboids() const;
  std::vector<Triangle<T> > const * getTriangles() const;

  std::vector<Sphere<T> > * getAndLockSpheres();
  std::vector<Cuboid<T> > * getAndLockCuboids();
  std::vector<Triangle<T> > * getAndLockTriangles();
  
 private:
  bool spheresLocked;
  bool cuboidsLocked;
  bool trianglesLocked;
  
  std::vector<Sphere<T> > spheres;
  std::vector<Cuboid<T> > cuboids;
  std::vector<Triangle<T> > triangles;

  void serializeMpiBroadcast(int root) const;
  void mpiBroadcastDeserialize(int root);

  void mpiBcastLargeArray(int root,
			  long long int arraySize,
			  double * array) const;
};

// ================================================================

/// Construct an empty model.
///
template <class T>
MixedModel<T>::MixedModel()
: spheresLocked(false),
  cuboidsLocked(false),
  trianglesLocked(false),
  spheres(),
  cuboids(),
  triangles() {

}

template <class T>
MixedModel<T>::~MixedModel() {

}

/// Add a sphere to the model.
///
template <class T>
void
MixedModel<T>::addSphere(Sphere<T> const & sphere) {
  assert(!spheresLocked);

  if (spheresLocked) {
    return;
  }
  
  if (sphere.getRadius() == 0) {
    std::cerr << "Warning: degenerate sphere" << std::endl
	      << sphere << std::endl;
  }
  
  spheres.push_back(sphere);
}

/// Add a cuboid to the model.
///
template <class T>
void
MixedModel<T>::addCuboid(Cuboid<T> const & cuboid) {
  assert(!cuboidsLocked);

  if (cuboidsLocked) {
    return;
  }
  
  if (cuboid.getMaxCoords() == cuboid.getMinCoords()) {
    std::cerr << "Warning: degenerate cuboid" << std::endl
	      << cuboid << std::endl;
  }
  
  cuboids.push_back(cuboid);
}


/// Add cuboids to the model representing voxels stored in a .fits.gz file.
///
template <class T>
bool
MixedModel<T>::loadVoxels(std::string const & inputFileName) {
  assert(!cuboidsLocked);

  if (cuboidsLocked) {
    return false;
  }
  
  Voxels<T> voxels;

  bool loadSuccess = voxels.loadFitsGz(inputFileName);

  if (!loadSuccess) {
    return false;
  }
  
  voxels.generateCuboids(&cuboids);
  
  return loadSuccess;
}

template <class T>
void
MixedModel<T>::mpiBroadcast(int root) {
#ifdef USE_MPI
  int mpiSize = 1, mpiRank = 0;

  MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

  if (mpiSize > 1) {
    if (mpiRank == root) {
      serializeMpiBroadcast(root);
    }
    else {
      mpiBroadcastDeserialize(root);
    }
  }
#endif
}

/// Broadcasts the geometric primitives stored in the model over MPI.
///
template <class T>
void
MixedModel<T>::serializeMpiBroadcast(int root) const {
#ifdef USE_MPI
  assert(!spheresLocked &&
	 !cuboidsLocked &&
	 !trianglesLocked);

  if (spheresLocked ||
      cuboidsLocked ||
      trianglesLocked) {
    return;
  }
  
  long long int sizeArray[3];

  sizeArray[0] = spheres.size();
  sizeArray[1] = cuboids.size();
  sizeArray[2] = triangles.size();
  
  MPI_Bcast(sizeArray, 3, MPI_LONG_LONG_INT, root, MPI_COMM_WORLD);

  long long int geometryArraySize =
    spheres.size() * 4 +
    cuboids.size() * 6 +
    triangles.size() * 9;
  
  double * geometryArray = new double[geometryArraySize];

  long long int geometryArrayIndex = 0;

  // spheres
  for (Sphere<T> const & sphere : spheres) {
    for (int dim = 0; dim < 3; ++ dim) {
      geometryArray[geometryArrayIndex++] = sphere.getCenter().get(dim);
    }
    
    geometryArray[geometryArrayIndex++] = sphere.getRadius();
  }

  // cuboids
  for (Cuboid<T> const & cuboid : cuboids) {
    for (int dim = 0; dim < 3; ++ dim) {
      geometryArray[geometryArrayIndex++] = cuboid.getMinCoords().get(dim);
    }

    for (int dim = 0; dim < 3; ++ dim) {
      geometryArray[geometryArrayIndex++] = cuboid.getMaxCoords().get(dim);
    }
  }

  // triangles
  for (Triangle<T> const & triangle : triangles) {
    for (int dim = 0; dim < 3; ++ dim) {
      geometryArray[geometryArrayIndex++] = triangle.getV1().get(dim);
    }

    for (int dim = 0; dim < 3; ++ dim) {
      geometryArray[geometryArrayIndex++] = triangle.getV2().get(dim);
    }

    for (int dim = 0; dim < 3; ++ dim) {
      geometryArray[geometryArrayIndex++] = triangle.getV3().get(dim);
    }
  }
  
  assert(geometryArrayIndex == geometryArraySize);

  mpiBcastLargeArray(root, geometryArraySize, geometryArray);

  delete [] geometryArray;
#endif
}

/// Receives a set of geometric primitives over MPI and stores them in the
/// model.
///
template <class T>
void
MixedModel<T>::mpiBroadcastDeserialize(int root) {
#ifdef USE_MPI
  assert(!spheresLocked &&
	 !cuboidsLocked &&
	 !trianglesLocked);

  if (spheresLocked ||
      cuboidsLocked ||
      trianglesLocked) {
    return;
  }
  
  long long int sizeArray[3];

  MPI_Bcast(&sizeArray, 3, MPI_LONG_LONG_INT, root, MPI_COMM_WORLD);

  long long int spheresSize   = sizeArray[0];
  long long int cuboidsSize   = sizeArray[1];
  long long int trianglesSize = sizeArray[2];
  
  long long int geometryArraySize =
    spheresSize * 4 +
    cuboidsSize * 6 +
    trianglesSize * 9;

  double * geometryArray = new double[geometryArraySize];

  mpiBcastLargeArray(root, geometryArraySize, geometryArray);

  spheres.clear();
  cuboids.clear();
  triangles.clear();

  spheres.reserve(spheresSize);
  cuboids.reserve(cuboidsSize);
  triangles.reserve(trianglesSize);
  
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

  // triangles
  for (long long int i = 0; i < trianglesSize; i++) {
    Vector3<T> v1;

    for (int dim = 0; dim < 3; ++ dim) {
      v1.set(dim, geometryArray[geometryArrayIndex++]);
    }

    Vector3<T> v2;

    for (int dim = 0; dim < 3; ++ dim) {
      v2.set(dim, geometryArray[geometryArrayIndex++]);
    }

    Vector3<T> v3;

    for (int dim = 0; dim < 3; ++ dim) {
      v3.set(dim, geometryArray[geometryArrayIndex++]);
    }

    triangles.emplace_back(v1, v2, v3);
  }

  assert(geometryArrayIndex == geometryArraySize);

  delete [] geometryArray;
#endif
}

/// Returns whether the model contains no geometric primitives.
///
template <class T>
bool
MixedModel<T>::isEmpty() const {
  return (spheres.empty() &&
	  cuboids.empty() &&
	  triangles.empty());
}

template <class T>
std::vector<Sphere<T> > const *
MixedModel<T>::getSpheres() const {
  return &spheres;
}

template <class T>
std::vector<Cuboid<T> > const *
MixedModel<T>::getCuboids() const {
  return &cuboids;
}

template <class T>
std::vector<Triangle<T> > const *
MixedModel<T>::getTriangles() const {
  return &triangles;
}

template <class T>
std::vector<Sphere<T> > *
MixedModel<T>::getAndLockSpheres() {
  assert(!spheresLocked);

  if (spheresLocked) {
    return nullptr;
  }

  spheresLocked = true;
  
  return &spheres;
}

template <class T>
std::vector<Cuboid<T> > *
MixedModel<T>::getAndLockCuboids() {
  assert(!cuboidsLocked);

  if (cuboidsLocked) {
    return nullptr;
  }

  cuboidsLocked = true;
  
  return &cuboids;
}

template <class T>
std::vector<Triangle<T> > *
MixedModel<T>::getAndLockTriangles() {
  assert(!trianglesLocked);

  if (trianglesLocked) {
    return nullptr;
  }

  trianglesLocked = true;
  
  return &triangles;
}

/// Uses MPI_Bcast to send or receive an array with size represented by
/// long long int instead of int
///
template <class T>
void
MixedModel<T>::mpiBcastLargeArray(int root,
				  long long int arraySize,
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
	      MPI_DOUBLE, root, MPI_COMM_WORLD);

    countBroadcastSoFar += countToBroadcast;
  }
#endif
}

}

#endif
