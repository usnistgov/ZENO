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
// Time-stamp: <2016-08-29 17:05:58 dcj>
//
// ================================================================

#include "ResultsZeno.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

/// Constructs the class to collect results with the given bounding sphere 
/// from the given number of threads,
/// and optionally save the hit point locations.
///
ResultsZeno::
ResultsZeno(Sphere<double> const & boundingSphere,
	    int numThreads,
	    bool saveHitPoints) 
  : boundingSphereRadius(boundingSphere.getRadius()),
    boundingSphereCenter(boundingSphere.getCenter()),
    numThreads(numThreads),
    saveHitPoints(saveHitPoints),
    numWalks(NULL),
    hitMissMean(NULL),
    hitMissM2(NULL),
    KPlus(NULL),
    KMinus(NULL),
    KPlusMean(NULL),
    KMinusMean(NULL),
    KPlusM2(NULL),
    KMinusM2(NULL),
    VPlus(NULL),
    VMinus(NULL),
    VPlusMean(NULL),
    VMinusMean(NULL),
    VPlusM2(NULL),
    VMinusM2(NULL),
    numWalksReduced(0),
    numHitsReduced(0),
    numHitsVarianceReduced(0),
    KPlusReduced(0, 0, 0),
    KMinusReduced(0, 0, 0),
    VPlusReduced(0, 0, 0, 
		 0, 0, 0, 
		 0, 0, 0),
    VMinusReduced(0, 0, 0, 
		  0, 0, 0, 
		  0, 0, 0),
    points(NULL),
    normals(NULL),
    charges(NULL),
    gatheredPoints(),
    gatheredNormals(),
    gatheredCharges(),
    reduced(true),
    hitPointsGathered(true) {

  numWalks = new double[numThreads];

  hitMissMean = new double[numThreads];

  hitMissM2 = new double[numThreads];

  KPlus  = new Vector3<double>[numThreads];
  KMinus = new Vector3<double>[numThreads];

  KPlusMean  = new Vector3<double>[numThreads];
  KMinusMean = new Vector3<double>[numThreads];  

  KPlusM2  = new Vector3<double>[numThreads];
  KMinusM2 = new Vector3<double>[numThreads];

  VPlus  = new Matrix3x3<double>[numThreads];
  VMinus = new Matrix3x3<double>[numThreads];

  VPlusMean  = new Matrix3x3<double>[numThreads];
  VMinusMean = new Matrix3x3<double>[numThreads];

  VPlusM2  = new Matrix3x3<double>[numThreads];
  VMinusM2 = new Matrix3x3<double>[numThreads];

  for (int threadNum = 0; threadNum < numThreads; threadNum++) {
    numWalks[threadNum] = 0;

    hitMissMean[threadNum] = 0;

    hitMissM2[threadNum] = 0;

    KPlus[threadNum].setXYZ(0, 0, 0);
    KMinus[threadNum].setXYZ(0, 0, 0);

    KPlusMean[threadNum].setXYZ(0, 0, 0);
    KMinusMean[threadNum].setXYZ(0, 0, 0);

    KPlusM2[threadNum].setXYZ(0, 0, 0);
    KMinusM2[threadNum].setXYZ(0, 0, 0);

    for (int component = 0; component < 3*3; component++) {
      VPlus[threadNum].set(component, 0);
      VMinus[threadNum].set(component, 0);

      VPlusMean[threadNum].set(component, 0);
      VMinusMean[threadNum].set(component, 0);

      VPlusM2[threadNum].set(component, 0);
      VMinusM2[threadNum].set(component, 0);
    }
  }

  points  = new std::vector<Vector3<double> >[numThreads];
  normals = new std::vector<Vector3<double> >[numThreads];
  charges = new std::vector<Vector3<char> >[numThreads];
}

ResultsZeno::
~ResultsZeno() {
  delete [] numWalks;

  delete [] hitMissMean;

  delete [] hitMissM2;

  delete [] KPlus;
  delete [] KMinus;

  delete [] KPlusMean;
  delete [] KMinusMean;

  delete [] KPlusM2;
  delete [] KMinusM2;

  delete [] VPlus;
  delete [] VMinus;

  delete [] VPlusMean;
  delete [] VMinusMean;

  delete [] VPlusM2;
  delete [] VMinusM2;

  delete [] points;
  delete [] normals;
  delete [] charges;
}

/// Record a miss from the given thread number.
///
void 
ResultsZeno::
recordMiss(int threadNum) {
  assert(threadNum >= 0 && threadNum < numThreads);

  reduced = false;

  double hitMissData = 0;

  Vector3<double> KPlusData(0, 0, 0);
  Vector3<double> KMinusData(0, 0, 0);

  Matrix3x3<double> VPlusData(0, 0, 0, 0, 0, 0, 0, 0, 0);
  Matrix3x3<double> VMinusData(0, 0, 0, 0, 0, 0, 0, 0, 0);

  numWalks[threadNum]++;

  updateVariance(threadNum,
		 hitMissData,
		 KPlusData, 
		 KMinusData,
		 VPlusData, 
		 VMinusData);
}

/// Perform a parallel reduction on the hit counts and other statistics and 
/// corresponding variances across threads and MPI nodes.
///
void 
ResultsZeno::
reduce() {
  if (reduced) {
    return;
  }

  numWalksReduced = 0;

  numHitsReduced = 0;

  numHitsVarianceReduced = 0;

  KPlusReduced.setXYZ(0, 0, 0);
  KMinusReduced.setXYZ(0, 0, 0);

  KPlusVarianceReduced.setXYZ(0, 0, 0);
  KMinusVarianceReduced.setXYZ(0, 0, 0);

  for (int component = 0; component < 3*3; component++) {
    VPlusReduced.set(component, 0);
    VMinusReduced.set(component, 0);

    VPlusVarianceReduced.set(component, 0);
    VMinusVarianceReduced.set(component, 0);
  }

  for (int threadNum = 0; threadNum < numThreads; threadNum++) {
    const double nn1 = (double)numWalks[threadNum] / (numWalks[threadNum] - 1);

    numWalksReduced += numWalks[threadNum];

    numHitsReduced += hitMissMean[threadNum] * numWalks[threadNum];

    numHitsVarianceReduced += hitMissM2[threadNum] * nn1;

    KPlusReduced  += KPlus[threadNum];
    KMinusReduced += KMinus[threadNum];

    KPlusVarianceReduced  += KPlusM2[threadNum] * nn1;
    KMinusVarianceReduced += KMinusM2[threadNum] * nn1;

    VPlusReduced  += VPlus[threadNum];
    VMinusReduced += VMinus[threadNum];

    VPlusVarianceReduced  += VPlusM2[threadNum] * nn1;
    VMinusVarianceReduced += VMinusM2[threadNum] * nn1;
  }

#ifdef USE_MPI
  const int mpiBufferSize = 51;

  double sendbuf[mpiBufferSize];

  int offset = 0;

  sendbuf[offset++] = numWalksReduced;

  sendbuf[offset++] = numHitsReduced;
  sendbuf[offset++] = numHitsVarianceReduced;

  for (int i = 0; i < 3; i++) {
    sendbuf[offset++] = KPlusReduced.get(i);
    sendbuf[offset++] = KPlusVarianceReduced.get(i);
  }

  for (int i = 0; i < 3; i++) {
    sendbuf[offset++] = KMinusReduced.get(i);
    sendbuf[offset++] = KMinusVarianceReduced.get(i);
  }

  for (int i = 0; i < 9; i++) {
    sendbuf[offset++] = VPlusReduced.get(i);
    sendbuf[offset++] = VPlusVarianceReduced.get(i);
  }

  for (int i = 0; i < 9; i++) {
    sendbuf[offset++] = VMinusReduced.get(i);
    sendbuf[offset++] = VMinusVarianceReduced.get(i);
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

  numWalksReduced = recvbuf[offset++];

  numHitsReduced         = recvbuf[offset++];
  numHitsVarianceReduced = recvbuf[offset++];

  for (int i = 0; i < 3; i++) {
    KPlusReduced.set(i, recvbuf[offset++]);
    KPlusVarianceReduced.set(i, recvbuf[offset++]);
  }

  for (int i = 0; i < 3; i++) {
    KMinusReduced.set(i, recvbuf[offset++]);
    KMinusVarianceReduced.set(i, recvbuf[offset++]);
  }

  for (int i = 0; i < 9; i++) {
    VPlusReduced.set(i, recvbuf[offset++]);
    VPlusVarianceReduced.set(i, recvbuf[offset++]);
  }

  for (int i = 0; i < 9; i++) {
    VMinusReduced.set(i, recvbuf[offset++]);
    VMinusVarianceReduced.set(i, recvbuf[offset++]);
  }
#endif

  reduced = true;
}

/// Gather the hit locations from all threads and MPI nodes.
///
void 
ResultsZeno::
gatherHitPoints() {
  if (hitPointsGathered) {
    return;
  }

  for (int threadNum = 0; threadNum < numThreads; threadNum++) {
    gatheredPoints.insert(gatheredPoints.end(), 
			 points[threadNum].begin(), 
			 points[threadNum].end());

    gatheredNormals.insert(gatheredNormals.end(),
			  normals[threadNum].begin(),
			  normals[threadNum].end());

    gatheredCharges.insert(gatheredCharges.end(),
			  charges[threadNum].begin(),
			  charges[threadNum].end());
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
  double * combinedNormalArray = new double[combinedArrayLength];
  char *   combinedChargeArray = new char[combinedArrayLength];

  int * combinedArrayOffsets = new int[mpiSize];

  combinedArrayOffsets[0] = 0;

  for (int i = 1; i < mpiSize; i++) {
    combinedArrayOffsets[i] = combinedArrayOffsets[i - 1] + arrayLengths[i - 1];
  }

  double * pointArray  = new double[arrayLengths[mpiRank]];
  double * normalArray = new double[arrayLengths[mpiRank]];
  char *   chargeArray = new char[arrayLengths[mpiRank]];

  for (int index = 0; index < arrayLengths[mpiRank]; index++) {
    pointArray[index]  = gatheredPoints[index / 3].get(index % 3);
    normalArray[index] = gatheredNormals[index / 3].get(index % 3);
    chargeArray[index] = gatheredCharges[index / 3].get(index % 3);
  }

  MPI_Allgatherv(pointArray, 
		 arrayLengths[mpiRank], MPI_DOUBLE,
		 combinedPointArray, 
		 arrayLengths, combinedArrayOffsets, MPI_DOUBLE,
		 MPI_COMM_WORLD);

  MPI_Allgatherv(normalArray, 
		 arrayLengths[mpiRank], MPI_DOUBLE,
		 combinedNormalArray, 
		 arrayLengths, combinedArrayOffsets, MPI_DOUBLE,
		 MPI_COMM_WORLD);

  MPI_Allgatherv(chargeArray,
		 arrayLengths[mpiRank], MPI_BYTE,
		 combinedChargeArray,
		 arrayLengths, combinedArrayOffsets, MPI_BYTE,
		 MPI_COMM_WORLD);

  gatheredPoints.clear();
  gatheredNormals.clear();
  gatheredCharges.clear();

  gatheredPoints.reserve(combinedArrayLength / 3);
  gatheredNormals.reserve(combinedArrayLength / 3);
  gatheredCharges.reserve(combinedArrayLength / 3);

  for (int index = 0; index < combinedArrayLength; index += 3) {
    gatheredPoints.emplace_back(combinedPointArray[index + 0],
				combinedPointArray[index + 1],
				combinedPointArray[index + 2]);

    gatheredNormals.emplace_back(combinedNormalArray[index + 0],
				 combinedNormalArray[index + 1],
				 combinedNormalArray[index + 2]);

    gatheredCharges.emplace_back(combinedChargeArray[index + 0],
				 combinedChargeArray[index + 1],
				 combinedChargeArray[index + 2]);
  }

  delete [] arrayLengths;

  delete [] combinedPointArray;
  delete [] combinedNormalArray;
  delete [] combinedChargeArray;

  delete [] combinedArrayOffsets;

  delete [] pointArray;
  delete [] normalArray;
  delete [] chargeArray;
#endif

  hitPointsGathered = true;
}

void 
ResultsZeno::
updateVariance(int threadNum,
	       double hitMissData,
	       Vector3<double> const & KPlusData, 
	       Vector3<double> const & KMinusData,
	       Matrix3x3<double> const & VPlusData, 
	       Matrix3x3<double> const & VMinusData) {

  updateItemVariance(hitMissData,
		     numWalks[threadNum],
		     &(hitMissMean[threadNum]),
		     &(hitMissM2[threadNum]));

  updateItemVariance(KPlusData,
		     numWalks[threadNum],
		     &(KPlusMean[threadNum]),
		     &(KPlusM2[threadNum]));

  updateItemVariance(KMinusData,
		     numWalks[threadNum],
		     &(KMinusMean[threadNum]),
		     &(KMinusM2[threadNum]));

  updateItemVariance(VPlusData,
		     numWalks[threadNum],
		     &(VPlusMean[threadNum]),
		     &(VPlusM2[threadNum]));

  updateItemVariance(VMinusData,
		     numWalks[threadNum],
		     &(VMinusMean[threadNum]),
		     &(VMinusM2[threadNum]));
}

double 
ResultsZeno::
getNumWalks() const {
  assert(reduced);

  return numWalksReduced;
}

Uncertain<double> 
ResultsZeno::
getNumHits() const {
  assert(reduced);

  return Uncertain<double>(numHitsReduced, numHitsVarianceReduced);
}

Vector3<Uncertain<double> > 
ResultsZeno::
getKPlus() const {
  assert(reduced);

  return Uncertain<double>::zip(KPlusReduced, KPlusVarianceReduced);
}

Vector3<Uncertain<double> > 
ResultsZeno::
getKMinus() const {
  assert(reduced);

  return Uncertain<double>::zip(KMinusReduced, KMinusVarianceReduced);
}

Matrix3x3<Uncertain<double> > 
ResultsZeno::
getVPlus() const {
  assert(reduced);

  return Uncertain<double>::zip(VPlusReduced, VPlusVarianceReduced);
}

Matrix3x3<Uncertain<double> > 
ResultsZeno::
getVMinus() const {
  assert(reduced);

  return Uncertain<double>::zip(VMinusReduced, VMinusVarianceReduced);
}

bool
ResultsZeno::
getSaveHitPoints() const {

  return saveHitPoints;
}

std::vector<Vector3<double> > const * 
ResultsZeno::
getPoints() const {
  assert(saveHitPoints);
  assert(hitPointsGathered);

  return &gatheredPoints;
}

std::vector<Vector3<double> > const * 
ResultsZeno:: 
getNormals() const {
  assert(saveHitPoints);
  assert(hitPointsGathered);

  return &gatheredNormals;
}

std::vector<Vector3<char> > const * 
ResultsZeno:: 
getCharges() const {
  assert(saveHitPoints);
  assert(hitPointsGathered);

  return &gatheredCharges;
}

// ================================================================

// Local Variables:
// time-stamp-line-limit: 30
// mode: c++
// End:
