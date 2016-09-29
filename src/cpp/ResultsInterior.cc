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
// Date:    Wed Apr 22 11:11:48 2015 EDT
//
// Time-stamp: <2016-08-29 16:20:36 dcj>
//
// ================================================================

#ifdef USE_MPI
#include <mpi.h>
#endif

#include <cassert>
#include <cstdlib>

#include "ResultsInterior.h"

// ================================================================

/// Constructs the class to collect results from the given number of threads,
/// and optionally save the hit point locations.
///
ResultsInterior::
ResultsInterior(int numThreads,
		bool saveHitPoints)
  : numThreads(numThreads),
    saveHitPoints(saveHitPoints),
    numSamples(NULL),
    hitMissMean(NULL),
    hitMissM2(NULL),
    hitPointsSqrMean(NULL),
    hitPointsSqrM2(NULL),
    hitPointsMean(NULL),
    hitPointsM2(NULL),
    numSamplesReduced(0),
    numHitsReduced(0),
    numHitsVarianceReduced(0),
    hitPointsSqrSumReduced(0, 0, 0, 
			   0, 0, 0, 
			   0, 0, 0),
    hitPointsSqrSumVarianceReduced(0, 0, 0, 
				   0, 0, 0, 
				   0, 0, 0),
    hitPointsSumReduced(0, 0, 0),
    hitPointsSumVarianceReduced(0, 0, 0),
    points(NULL),
    gatheredPoints(),
    reduced(true),
    hitPointsGathered(true) {

  numSamples = new double[numThreads];

  hitMissMean = new double[numThreads];
  hitMissM2   = new double[numThreads];

  hitPointsSqrMean = new Matrix3x3<double>[numThreads];
  hitPointsSqrM2   = new Matrix3x3<double>[numThreads];

  hitPointsMean = new Vector3<double>[numThreads];
  hitPointsM2   = new Vector3<double>[numThreads];

  for (int threadNum = 0; threadNum < numThreads; threadNum++) {
    numSamples[threadNum] = 0;

    hitMissMean[threadNum] = 0;
    hitMissM2[threadNum]   = 0;

    for (int component = 0; component < 3*3; component++) {
      hitPointsSqrMean[threadNum].set(component, 0);
      hitPointsSqrM2[threadNum].set(component, 0);
    }

    hitPointsMean[threadNum].setXYZ(0, 0, 0);
    hitPointsM2[threadNum].setXYZ(0, 0, 0);
  }

  points = new std::vector<Vector3<double> >[numThreads];
}

ResultsInterior::
~ResultsInterior() {
  delete [] numSamples;

  delete [] hitMissMean;
  delete [] hitMissM2;

  delete [] hitPointsSqrMean;
  delete [] hitPointsSqrM2;

  delete [] hitPointsMean;
  delete [] hitPointsM2;
}

/// Record a hit from the given thread number at the given location.
///
void 
ResultsInterior::
recordHit(int threadNum,
	  Vector3<double> const & point) {

  assert(threadNum >= 0 && threadNum < numThreads);

  reduced = false;

  numSamples[threadNum] ++;

  double hitMissData = 1;

  Matrix3x3<double> hitPointsSqrData;

  for (int row = 0; row < 3; ++row) {
    for (int col = 0; col < 3; ++col) {
      double element = point.get(row) * point.get(col);

      hitPointsSqrData.set(row, col, element);
    }
  }

  Vector3<double> hitPointsData(point);

  updateVariance(threadNum,
		 hitMissData,
		 hitPointsSqrData,
		 hitPointsData);

  if (saveHitPoints) {
    hitPointsGathered = false;

    points[threadNum].push_back(point);
  }
}

/// Record a miss from the given thread number.
///
void 
ResultsInterior::
recordMiss(int threadNum) {
  assert(threadNum >= 0 && threadNum < numThreads);

  reduced = false;

  numSamples[threadNum] ++;

  double hitMissData = 0;

  Matrix3x3<double> hitPointsSqrData(0, 0, 0,
				     0, 0, 0,
				     0, 0, 0);

  Vector3<double> hitPointsData(0, 0, 0);

  updateVariance(threadNum, 
		 hitMissData,
		 hitPointsSqrData,
		 hitPointsData);
}

/// Perform a parallel reduction on the hit counts and locations and 
/// corresponding variances across threads and MPI nodes.
///
void 
ResultsInterior::
reduce() {
  if (reduced) {
    return;
  }

  numSamplesReduced = 0;

  numHitsReduced         = 0;
  numHitsVarianceReduced = 0;

  for (int component = 0; component < 3*3; component++) {
    hitPointsSqrSumReduced.set(component, 0);
    hitPointsSqrSumVarianceReduced.set(component, 0);
  }

  hitPointsSumReduced.setXYZ(0, 0, 0);
  hitPointsSumVarianceReduced.setXYZ(0, 0, 0);

  for (int threadNum = 0; threadNum < numThreads; threadNum++) {

    numSamplesReduced += numSamples[threadNum];

    reduceItem(hitMissMean[threadNum],
	       hitMissM2[threadNum],
	       numSamples[threadNum],
	       &numHitsReduced,
	       &numHitsVarianceReduced);

    reduceItem(hitPointsSqrMean[threadNum],
	       hitPointsSqrM2[threadNum],
	       numSamples[threadNum],
	       &hitPointsSqrSumReduced,
	       &hitPointsSqrSumVarianceReduced);

    reduceItem(hitPointsMean[threadNum],
	       hitPointsM2[threadNum],
	       numSamples[threadNum],
	       &hitPointsSumReduced,
	       &hitPointsSumVarianceReduced);
  }

#ifdef USE_MPI
  const int mpiBufferSize = 27; 

  double sendbuf[mpiBufferSize];

  int offset = 0;

  sendbuf[offset++] = numSamplesReduced;

  sendbuf[offset++] = numHitsReduced;
  sendbuf[offset++] = numHitsVarianceReduced;

  for (int i = 0; i < 9; i++) {
    sendbuf[offset++] = hitPointsSqrSumReduced.get(i);
    sendbuf[offset++] = hitPointsSqrSumVarianceReduced.get(i);
  }

  for (int i = 0; i < 3; i++) {
    sendbuf[offset++] = hitPointsSumReduced.get(i);
    sendbuf[offset++] = hitPointsSumVarianceReduced.get(i);
  }

  double recvbuf[mpiBufferSize];

  for (int i = 0; i < mpiBufferSize; i++) {
    recvbuf[i] = 0;
  }

  // MPI_Reduce(sendbuf, recvbuf, mpiBufferSize, MPI_DOUBLE,
  // 	     MPI_SUM, 0, MPI_COMM_WORLD);

  MPI_Allreduce(sendbuf, recvbuf, mpiBufferSize, MPI_DOUBLE,
		MPI_SUM, MPI_COMM_WORLD);

  offset = 0;

  numSamplesReduced = recvbuf[offset++];

  numHitsReduced         = recvbuf[offset++];
  numHitsVarianceReduced = recvbuf[offset++];

  for (int i = 0; i < 9; i++) {
    hitPointsSqrSumReduced.set(i, recvbuf[offset++]);
    hitPointsSqrSumVarianceReduced.set(i, recvbuf[offset++]);
  }

  for (int i = 0; i < 3; i++) {
    hitPointsSumReduced.set(i, recvbuf[offset++]);
    hitPointsSumVarianceReduced.set(i, recvbuf[offset++]);
  }
#endif

  reduced = true;
}

/// Gather the hit locations from all threads and MPI nodes.
///
void 
ResultsInterior::
gatherHitPoints() {
  if (hitPointsGathered) {
    return;
  }

  for (int threadNum = 0; threadNum < numThreads; threadNum++) {
    gatheredPoints.insert(gatheredPoints.end(), 
			 points[threadNum].begin(), 
			 points[threadNum].end());
  }

#ifdef USE_MPI
  int mpiSize = 0;
  int mpiRank = 0;

  MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

  int * arrayLengths = new int[mpiSize];

  int arrayLength = gatheredPoints.size() * 3;

  MPI_Allgather(&arrayLength, 1, MPI_INT,
		arrayLengths, 1, MPI_INT,
		MPI_COMM_WORLD);

  int combinedArrayLength = 0;

  for (int i = 0; i < mpiSize; i++) {
    combinedArrayLength += arrayLengths[i];
  }

  double * combinedPointArray  = new double[combinedArrayLength];

  int * combinedArrayOffsets = new int[mpiSize];

  combinedArrayOffsets[0] = 0;

  for (int i = 1; i < mpiSize; i++) {
    combinedArrayOffsets[i] = combinedArrayOffsets[i - 1] + arrayLengths[i - 1];
  }

  double * pointArray  = new double[arrayLengths[mpiRank]];

  for (int index = 0; index < arrayLengths[mpiRank]; index++) {
    pointArray[index] = gatheredPoints[index / 3].get(index % 3);
  }

  MPI_Allgatherv(pointArray, 
		 arrayLengths[mpiRank], MPI_DOUBLE,
		 combinedPointArray, 
		 arrayLengths, combinedArrayOffsets, MPI_DOUBLE,
		 MPI_COMM_WORLD);

  gatheredPoints.clear();

  gatheredPoints.reserve(combinedArrayLength / 3);
  
  for (int index = 0; index < combinedArrayLength; index += 3) {
    gatheredPoints.emplace_back(combinedPointArray[index + 0],
				combinedPointArray[index + 1],
				combinedPointArray[index + 2]);
  }

  delete [] arrayLengths;

  delete [] combinedPointArray;

  delete [] combinedArrayOffsets;

  delete [] pointArray;
#endif

  hitPointsGathered = true;
}

Uncertain<double> 
ResultsInterior::
getNumHits() const {
  assert(reduced);

  return Uncertain<double>(numHitsReduced, numHitsVarianceReduced);
}

Matrix3x3<Uncertain<double> > 
ResultsInterior::
getHitPointsSqrSum() const {
  assert(reduced);

  return Uncertain<double>::zip(hitPointsSqrSumReduced, 
				hitPointsSqrSumVarianceReduced);
}

Vector3<Uncertain<double> > 
ResultsInterior::
getHitPointsSum() const {
  assert(reduced);

  return Uncertain<double>::zip(hitPointsSumReduced, 
				hitPointsSumVarianceReduced);
}

double
ResultsInterior::
getNumSamples() const {
  assert(reduced);

  return numSamplesReduced;
}

bool
ResultsInterior::
getSaveHitPoints() const {

  return saveHitPoints;
}

std::vector<Vector3<double> > const * 
ResultsInterior::
getPoints() const {
  assert(saveHitPoints);
  assert(hitPointsGathered);

  return &gatheredPoints;
}

void 
ResultsInterior::
updateVariance(int threadNum,
	       double hitMissData,
	       Matrix3x3<double> const & hitPointsSqrData,
	       Vector3<double> const & hitPointsData) {

  updateItemVariance(hitMissData,
		     numSamples[threadNum],
		     &(hitMissMean[threadNum]),
		     &(hitMissM2[threadNum]));

  updateItemVariance(hitPointsSqrData,
		     numSamples[threadNum],
		     &(hitPointsSqrMean[threadNum]),
		     &(hitPointsSqrM2[threadNum]));

  updateItemVariance(hitPointsData,
		     numSamples[threadNum],
		     &(hitPointsMean[threadNum]),
		     &(hitPointsM2[threadNum]));
}

// ================================================================

// Local Variables:
// time-stamp-line-limit: 30
// mode: c++
// End:
