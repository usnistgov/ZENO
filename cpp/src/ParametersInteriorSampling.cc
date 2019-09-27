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

#include "ParametersInteriorSampling.h"

// ================================================================

using namespace zeno;

ParametersInteriorSampling::ParametersInteriorSampling() 
  : numThreads(),
    seed(),
    totalNumSamples(),
    totalNumSamplesWasSet(false),
    maxErrorVolume(),
    maxErrorVolumeWasSet(false),
    maxRunTime(),
    maxRunTimeWasSet(false),
    minTotalNumSamples(),
    saveInteriorPoints(false),
    launchCenter(),
    launchCenterWasSet(false),
    launchRadius(),
    launchRadiusWasSet(false) {

}

ParametersInteriorSampling::~ParametersInteriorSampling() {

}

void
ParametersInteriorSampling::setNumThreads(int numThreads) {
  this->numThreads = numThreads;
}

int 
ParametersInteriorSampling::getNumThreads() const {
  return numThreads;
}

void
ParametersInteriorSampling::setSeed(int seed) {
  this->seed = seed;
}

int 
ParametersInteriorSampling::getSeed() const {
  return seed;
}

void
ParametersInteriorSampling::setTotalNumSamples(long long totalNumSamples) {
  this->totalNumSamples = totalNumSamples;

  totalNumSamplesWasSet = true;
}

long long 
ParametersInteriorSampling::getTotalNumSamples() const {
  return totalNumSamples;
}

bool 
ParametersInteriorSampling::getTotalNumSamplesWasSet() const {
  return totalNumSamplesWasSet;
}

void
ParametersInteriorSampling::setMaxRunTime(double maxRunTime) {
  this->maxRunTime = maxRunTime;

  maxRunTimeWasSet = true;
}

double
ParametersInteriorSampling::getMaxRunTime() const {
  return maxRunTime;
}

bool
ParametersInteriorSampling::getMaxRunTimeWasSet() const {
  return maxRunTimeWasSet;
}

void
ParametersInteriorSampling::setMaxErrorVolume(double maxErrorVolume) {
  this->maxErrorVolume = maxErrorVolume;

  maxErrorVolumeWasSet = true;
}

double
ParametersInteriorSampling::getMaxErrorVolume() const {
  return maxErrorVolume;
}

bool
ParametersInteriorSampling::getMaxErrorVolumeWasSet() const {
  return maxErrorVolumeWasSet;
}

void
ParametersInteriorSampling::setMinTotalNumSamples(long long minTotalNumSamples) {
  this->minTotalNumSamples = minTotalNumSamples;
}

long long 
ParametersInteriorSampling::getMinTotalNumSamples() const {
  return minTotalNumSamples;
}

void
ParametersInteriorSampling::setSaveInteriorPoints(bool saveInteriorPoints) {
  this->saveInteriorPoints = saveInteriorPoints;
}

bool
ParametersInteriorSampling::getSaveInteriorPoints() const {
  return saveInteriorPoints;
}

void 
ParametersInteriorSampling::setLaunchCenter(Vector3<double> launchCenter) {
  this->launchCenter = launchCenter;

  launchCenterWasSet = true;
}

Vector3<double>
ParametersInteriorSampling::getLaunchCenter() const {
  return launchCenter;
}

bool 
ParametersInteriorSampling::getLaunchCenterWasSet() const {
  return launchCenterWasSet;
}

void 
ParametersInteriorSampling::setLaunchRadius(double launchRadius) {
  this->launchRadius = launchRadius;

  launchRadiusWasSet = true;
}

double
ParametersInteriorSampling::getLaunchRadius() const {
  return launchRadius;
}

bool 
ParametersInteriorSampling::getLaunchRadiusWasSet() const {
  return launchRadiusWasSet;
}

Sphere<double>
ParametersInteriorSampling::getLaunchSphere() const {
  return Sphere<double>(launchCenter,
			launchRadius);
}

void 
ParametersInteriorSampling::mpiBroadcast(int root) {
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
ParametersInteriorSampling::serializeMpiBroadcast(int root) const {
#ifdef USE_MPI
  const int numDoublesToSend = 10;
  
  double doublesArray[numDoublesToSend];

  doublesArray[0] = getMaxErrorVolume();
  doublesArray[1] = (double)getMaxErrorVolumeWasSet();
  doublesArray[2] = getMaxRunTime();
  doublesArray[3] = (double)getMaxRunTimeWasSet();
  doublesArray[4] = getLaunchCenter().getX();
  doublesArray[5] = getLaunchCenter().getY();
  doublesArray[6] = getLaunchCenter().getZ();
  doublesArray[7] = (double)getLaunchCenterWasSet();
  doublesArray[8] = getLaunchRadius();
  doublesArray[9] = (double)getLaunchRadiusWasSet();

  MPI_Bcast(doublesArray, numDoublesToSend, MPI_DOUBLE,
	    root, MPI_COMM_WORLD);
  
  const int numLongLongsToSend = 6;

  long long longLongsArray[numLongLongsToSend];

  longLongsArray[0] = (long long)getNumThreads();
  longLongsArray[1] = (long long)getSeed();
  longLongsArray[2] = getTotalNumSamples();
  longLongsArray[3] = (long long)getTotalNumSamplesWasSet();
  longLongsArray[4] = getMinTotalNumSamples();
  longLongsArray[5] = (long long)getSaveInteriorPoints();

  MPI_Bcast(longLongsArray, numLongLongsToSend, MPI_LONG_LONG,
	    root, MPI_COMM_WORLD);
#endif
}

/// Receives the parameters that cannot be set in the bod file over MPI.
///
void 
ParametersInteriorSampling::mpiBroadcastDeserialize(int root) {
#ifdef USE_MPI
  const int numDoublesToReceive = 10;

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
    setTotalNumSamples(longLongsArray[2]);
  }

  if ((bool)doublesArray[1]) {
    setMaxErrorVolume(doublesArray[0]);
  }

  if ((bool)doublesArray[3]) {
    setMaxRunTime(doublesArray[2]);
  }

  setMinTotalNumSamples(longLongsArray[4]);

  setSaveInteriorPoints((bool)longLongsArray[5]);

  if ((bool)doublesArray[7]) {
    setLaunchCenter(Vector3<double>(doublesArray[4],
				    doublesArray[5],
				    doublesArray[6]));
  }

  if ((bool)doublesArray[9]) {
    setLaunchRadius(doublesArray[8]);
  }
#endif
}
