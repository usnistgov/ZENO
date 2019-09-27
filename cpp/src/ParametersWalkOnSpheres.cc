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
// Created: 2019-07-03
// 
// ================================================================

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "ParametersWalkOnSpheres.h"

// ================================================================

using namespace zeno;

ParametersWalkOnSpheres::ParametersWalkOnSpheres() 
  : numThreads(),
    seed(),
    totalNumWalks(),
    totalNumWalksWasSet(false),
    maxErrorCapacitance(),
    maxErrorCapacitanceWasSet(false),
    maxErrorPolarizability(),
    maxErrorPolarizabilityWasSet(false),
    maxRunTime(),
    maxRunTimeWasSet(false),
    minTotalNumWalks(),
    saveSurfacePoints(false),
    skinThickness(),
    skinThicknessWasSet(false),
    launchCenter(),
    launchCenterWasSet(false),
    launchRadius(),
    launchRadiusWasSet(false) {

}

ParametersWalkOnSpheres::~ParametersWalkOnSpheres() {

}

void
ParametersWalkOnSpheres::setNumThreads(int numThreads) {
  this->numThreads = numThreads;
}

int 
ParametersWalkOnSpheres::getNumThreads() const {
  return numThreads;
}

void
ParametersWalkOnSpheres::setSeed(int seed) {
  this->seed = seed;
}

int 
ParametersWalkOnSpheres::getSeed() const {
  return seed;
}

void
ParametersWalkOnSpheres::setTotalNumWalks(long long totalNumWalks) {
  this->totalNumWalks = totalNumWalks;

  totalNumWalksWasSet = true;
}

long long 
ParametersWalkOnSpheres::getTotalNumWalks() const {
  return totalNumWalks;
}

bool
ParametersWalkOnSpheres::getTotalNumWalksWasSet() const {
  return totalNumWalksWasSet;
}

void
ParametersWalkOnSpheres::setMaxErrorCapacitance(double maxErrorCapacitance) {
  this->maxErrorCapacitance = maxErrorCapacitance;

  maxErrorCapacitanceWasSet = true;
}

double
ParametersWalkOnSpheres::getMaxErrorCapacitance() const {
  return maxErrorCapacitance;
}

bool
ParametersWalkOnSpheres::getMaxErrorCapacitanceWasSet() const {
  return maxErrorCapacitanceWasSet;
}

void
ParametersWalkOnSpheres::setMaxErrorPolarizability
(double maxErrorPolarizability) {
  this->maxErrorPolarizability = maxErrorPolarizability;

  maxErrorPolarizabilityWasSet = true;
}

double
ParametersWalkOnSpheres::getMaxErrorPolarizability() const {
  return maxErrorPolarizability;
}

bool
ParametersWalkOnSpheres::getMaxErrorPolarizabilityWasSet() const {
  return maxErrorPolarizabilityWasSet;
}

void
ParametersWalkOnSpheres::setMaxRunTime(double maxRunTime) {
  this->maxRunTime = maxRunTime;

  maxRunTimeWasSet = true;
}

double
ParametersWalkOnSpheres::getMaxRunTime() const {
  return maxRunTime;
}

bool
ParametersWalkOnSpheres::getMaxRunTimeWasSet() const {
  return maxRunTimeWasSet;
}

void
ParametersWalkOnSpheres::setMinTotalNumWalks(long long minTotalNumWalks) {
  this->minTotalNumWalks = minTotalNumWalks;
}

long long 
ParametersWalkOnSpheres::getMinTotalNumWalks() const {
  return minTotalNumWalks;
}

void
ParametersWalkOnSpheres::setSaveSurfacePoints(bool saveSurfacePoints) {
  this->saveSurfacePoints = saveSurfacePoints;
}

bool
ParametersWalkOnSpheres::getSaveSurfacePoints() const {
  return saveSurfacePoints;
}

void 
ParametersWalkOnSpheres::setSkinThickness(double skinThickness) {
  this->skinThickness = skinThickness;

  skinThicknessWasSet = true;
}

double
ParametersWalkOnSpheres::getSkinThickness() const {
  return skinThickness;
}

bool
ParametersWalkOnSpheres::getSkinThicknessWasSet() const {
  return skinThicknessWasSet;
}

void 
ParametersWalkOnSpheres::setLaunchCenter(Vector3<double> launchCenter) {
  this->launchCenter = launchCenter;

  launchCenterWasSet = true;
}

Vector3<double>
ParametersWalkOnSpheres::getLaunchCenter() const {
  return launchCenter;
}

bool 
ParametersWalkOnSpheres::getLaunchCenterWasSet() const {
  return launchCenterWasSet;
}

void 
ParametersWalkOnSpheres::setLaunchRadius(double launchRadius) {
  this->launchRadius = launchRadius;

  launchRadiusWasSet = true;
}

double
ParametersWalkOnSpheres::getLaunchRadius() const {
  return launchRadius;
}

bool 
ParametersWalkOnSpheres::getLaunchRadiusWasSet() const {
  return launchRadiusWasSet;
}

Sphere<double>
ParametersWalkOnSpheres::getLaunchSphere() const {
  return Sphere<double>(launchCenter,
			launchRadius);
}

void 
ParametersWalkOnSpheres::mpiBroadcast(int root) {
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

/// Broadcasts the parameters that cannot be set in the bod file over MPI.
///
void 
ParametersWalkOnSpheres::serializeMpiBroadcast(int root) const {
#ifdef USE_MPI
  const int numDoublesToSend = 14;
  
  double doublesArray[numDoublesToSend];

  doublesArray[0]  = getMaxErrorCapacitance();
  doublesArray[1]  = (double)getMaxErrorCapacitanceWasSet();
  doublesArray[2]  = getMaxErrorPolarizability();
  doublesArray[3]  = (double)getMaxErrorPolarizabilityWasSet();
  doublesArray[4]  = getMaxRunTime();
  doublesArray[5]  = (double)getMaxRunTimeWasSet();
  doublesArray[6]  = getSkinThickness();
  doublesArray[7]  = (double)getSkinThicknessWasSet();
  doublesArray[8]  = getLaunchCenter().getX();
  doublesArray[9]  = getLaunchCenter().getY();
  doublesArray[10] = getLaunchCenter().getZ();
  doublesArray[11] = (double)getLaunchCenterWasSet();
  doublesArray[12] = getLaunchRadius();
  doublesArray[13] = (double)getLaunchRadiusWasSet();

  MPI_Bcast(doublesArray, numDoublesToSend, MPI_DOUBLE,
	    root, MPI_COMM_WORLD);
  
  const int numLongLongsToSend = 6;

  long long longLongsArray[numLongLongsToSend];

  longLongsArray[0] = (long long)getNumThreads();
  longLongsArray[1] = (long long)getSeed();
  longLongsArray[2] = getTotalNumWalks();
  longLongsArray[3] = (long long)getTotalNumWalksWasSet();
  longLongsArray[4] = getMinTotalNumWalks();
  longLongsArray[5] = (long long)getSaveSurfacePoints();

  MPI_Bcast(longLongsArray, numLongLongsToSend, MPI_LONG_LONG,
	    root, MPI_COMM_WORLD);
#endif
}

/// Receives the parameters that cannot be set in the bod file over MPI.
///
void 
ParametersWalkOnSpheres::mpiBroadcastDeserialize(int root) {
#ifdef USE_MPI
  const int numDoublesToReceive = 14;

  double doublesArray[numDoublesToReceive];

  MPI_Bcast(doublesArray, numDoublesToReceive, MPI_DOUBLE,
	    root, MPI_COMM_WORLD);

  const int numLongLongsToReceive = 6;

  long long longLongsArray[numLongLongsToReceive];

  MPI_Bcast(longLongsArray, numLongLongsToReceive, MPI_LONG_LONG,
	    root, MPI_COMM_WORLD);

  setNumThreads((int)longLongsArray[0]);

  setSeed((int)longLongsArray[1]);

  if ((bool)longLongsArray[3]) {
    setTotalNumWalks(longLongsArray[2]);
  }

  if ((bool)doublesArray[1]) {
    setMaxErrorCapacitance(doublesArray[0]);
  }

  if ((bool)doublesArray[3]) {
    setMaxErrorPolarizability(doublesArray[2]);
  }

  if ((bool)doublesArray[5]) {
    setMaxRunTime(doublesArray[4]);
  }

  setMinTotalNumWalks(longLongsArray[4]);

  setSaveSurfacePoints((bool)longLongsArray[5]);

  if ((bool)doublesArray[7]) {
    setSkinThickness(doublesArray[6]);
  }

  if ((bool)doublesArray[11]) {
    setLaunchCenter(Vector3<double>(doublesArray[8],
				    doublesArray[9],
				    doublesArray[10]));
  }

  if ((bool)doublesArray[13]) {
    setLaunchRadius(doublesArray[12]);
  }
#endif
}
