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
// Time-stamp: <2016-09-19 16:04:34 dcj>
//
// ================================================================

#ifndef SPHERES_H
#define SPHERES_H

#ifdef USE_MPI
#include <mpi.h>
#endif

#include <vector>

#include "Sphere.h"

// ================================================================

/// Represents a set of spheres.
///
template <class T>
class Spheres {
public:
  Spheres();
  ~Spheres();

  void add(Sphere<T> const & sphere);

  std::vector<Sphere<T> > const & getVector() const;

  bool isEmpty() const;

  void mpiSend() const;
  void mpiReceive();

private:
  std::vector<Sphere<T> > spheres;
};

/// Construct an empty set of spheres.
///
template <class T>
Spheres<T>::Spheres() 
  : spheres() {

}

template <class T>
Spheres<T>::~Spheres() {

}

/// Add a sphere to the set.
///
template <class T>
void 
Spheres<T>::add(Sphere<T> const & sphere) {

  spheres.push_back(sphere);
}

/// Returns the set of spheres as a vector.
///
template <class T>
std::vector<Sphere<T> > const &
Spheres<T>::getVector() const {

  return spheres;
}

template <class T>
bool 
Spheres<T>::isEmpty() const {

  return (spheres.size() == 0);
}

/// Broadcasts the set of spheres over MPI.
///
template <class T>
void
Spheres<T>::mpiSend() const {
#ifdef USE_MPI
  int numSpheres = spheres.size();

  MPI_Bcast(&numSpheres, 1, MPI_INT, 0, MPI_COMM_WORLD);

  double * sphereArray = new double[numSpheres*4];

  for (int i = 0; i < numSpheres; i++) {
    sphereArray[i*4 + 0] = spheres.at(i).getCenter().get(0);
    sphereArray[i*4 + 1] = spheres.at(i).getCenter().get(1);
    sphereArray[i*4 + 2] = spheres.at(i).getCenter().get(2);
    sphereArray[i*4 + 3] = spheres.at(i).getRadius();
  }

  MPI_Bcast(sphereArray, numSpheres*4, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  delete [] sphereArray;
#endif
}

/// Receives a set of spheres over MPI and adds it to the current set.
///
template <class T>
void
Spheres<T>::mpiReceive() {
#ifdef USE_MPI
  int numSpheres = 0;

  MPI_Bcast(&numSpheres, 1, MPI_INT, 0, MPI_COMM_WORLD);

  double * sphereArray = new double[numSpheres * 4];

  MPI_Bcast(sphereArray, numSpheres*4, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  for (int i = 0; i < numSpheres; i++) {
    spheres.emplace_back(Vector3<double>(sphereArray[i*4 + 0],
					 sphereArray[i*4 + 1],
					 sphereArray[i*4 + 2]),
			 sphereArray[i*4 + 3]);
  }

  delete [] sphereArray;
#endif
}

// ================================================================

#endif

// ================================================================

// Local Variables:
// time-stamp-line-limit: 30
// End:
