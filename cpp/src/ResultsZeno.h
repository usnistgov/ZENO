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
// Created: Wed Apr 22 11:11:48 2015 EDT
//
// ================================================================

#ifndef RESULTS_ZENO_H_
#define RESULTS_ZENO_H_

// ================================================================

#include <iostream>
#include <algorithm>
#include <cmath>
#include <cassert>

#include "Geometry/Vector3.h"
#include "Geometry/Matrix3x3.h"
#include "Geometry/Sphere.h"

#include "Uncertain.h"

// ================================================================

namespace zeno {

/// Collects results from the Walk-on-Spheres computation.
///
class ResultsZeno {
public:
  ResultsZeno(Sphere<double> const & boundingSphere,
	      int numThreads,
	      bool saveHitPoints);

  ~ResultsZeno();

  template <class RandomNumberGenerator>
  void recordHit(int threadNum,
		 Vector3<double> const & startPoint,
		 Vector3<double> const & endPoint,
		 RandomNumberGenerator * randomNumberGenerator);

  void recordMiss(int threadNum);

  void reduce();

  void gatherHitPoints();

  double getNumWalks() const;

  Uncertain<double> getNumHits() const;

  Vector3<Uncertain<double> > getKPlus() const;
  Vector3<Uncertain<double> > getKMinus() const;

  Matrix3x3<Uncertain<double> > getVPlus() const;
  Matrix3x3<Uncertain<double> > getVMinus() const;

  bool getSaveHitPoints() const;

  Sphere<double> getBoundingSphere() const;

  std::vector<Vector3<double> > const * getPoints() const;
  std::vector<Vector3<char> > const * getCharges() const;

private:
  void updateVariance(int threadNum,
		      double hitMissData,
		      Vector3<double> const & KPlusData, 
		      Vector3<double> const & KMinusData,
		      Matrix3x3<double> const & VPlusData, 
		      Matrix3x3<double> const & VMinusData);

  template <class T>
  void updateItemVariance(T const & data,
			  double num,
			  T * mean,
			  T * M2);

  Sphere<double> const boundingSphere;

  int const numThreads;

  bool saveHitPoints;

  double * numWalks;

  double * hitMissMean;
  double * hitMissM2;

  Vector3<double> * KPlus;
  Vector3<double> * KMinus;

  Vector3<double> * KPlusMean;
  Vector3<double> * KMinusMean;

  Vector3<double> * KPlusM2;
  Vector3<double> * KMinusM2;

  Matrix3x3<double> * VPlus;
  Matrix3x3<double> * VMinus;

  Matrix3x3<double> * VPlusMean;
  Matrix3x3<double> * VMinusMean;

  Matrix3x3<double> * VPlusM2;
  Matrix3x3<double> * VMinusM2;

  double numWalksReduced;

  double numHitsReduced;
  double numHitsVarianceReduced;

  Vector3<double> KPlusReduced;
  Vector3<double> KMinusReduced;

  Vector3<double> KPlusVarianceReduced;
  Vector3<double> KMinusVarianceReduced;

  Matrix3x3<double> VPlusReduced;
  Matrix3x3<double> VMinusReduced;

  Matrix3x3<double> VPlusVarianceReduced;
  Matrix3x3<double> VMinusVarianceReduced;

  std::vector<Vector3<double> > * points;
  std::vector<Vector3<char> > * charges;

  std::vector<Vector3<double> > gatheredPoints;
  std::vector<Vector3<char> > gatheredCharges;

  bool reduced;
  bool hitPointsGathered;
};

/// Record a hit from the given thread number of a walk from the given start
/// point to the given end point with the given end point surface normal.
/// Walker charges are assigned using the given random number generator.
///
template <class RandomNumberGenerator>
void 
ResultsZeno::
recordHit(int threadNum,
	  Vector3<double> const & startPoint,
	  Vector3<double> const & endPoint,
	  RandomNumberGenerator * randomNumberGenerator) {

  assert(threadNum >= 0 && threadNum < numThreads);

  reduced = false;

  double const boundingSphereRadius = boundingSphere.getRadius();
  Vector3<double> const boundingSphereCenter = boundingSphere.getCenter();

  Vector3<double> normalizedStartPoint = startPoint - boundingSphereCenter;
  Vector3<double> normalizedEndPoint   = endPoint - boundingSphereCenter;

  Vector3<char> walkCharges;

  double hitMissData = 1;

  Vector3<double> KPlusData(0, 0, 0);
  Vector3<double> KMinusData(0, 0, 0);

  Matrix3x3<double> VPlusData(0, 0, 0, 0, 0, 0, 0, 0, 0);
  Matrix3x3<double> VMinusData(0, 0, 0, 0, 0, 0, 0, 0, 0);

  numWalks[threadNum]++;

  for (int dim = 0; dim < 3; dim++) {
    double probability = 
      0.5 + normalizedStartPoint.get(dim)/(2*boundingSphereRadius);

    if (probability > randomNumberGenerator->getRandIn01()) {
      //c[dim] == +1

      walkCharges.set(dim, '+');

      KPlusData.add(dim, 1);
      VPlusData.addRow(dim, normalizedEndPoint);
    }
    else {
      //c[dim] == -1

      walkCharges.set(dim, '-');

      KMinusData.add(dim, 1);
      VMinusData.addRow(dim, normalizedEndPoint);
    }
  }

  KPlus[threadNum]  += KPlusData;
  KMinus[threadNum] += KMinusData;

  VPlus[threadNum]  += VPlusData;
  VMinus[threadNum] += VMinusData;

  updateVariance(threadNum,
		 hitMissData,
		 KPlusData, 
		 KMinusData,
		 VPlusData, 
		 VMinusData);

  if (saveHitPoints) {
    hitPointsGathered = false;

    points[threadNum].push_back(endPoint);
    charges[threadNum].push_back(walkCharges);
  }
}

template <class T>
void
ResultsZeno::
updateItemVariance(T const & data,
		   double num,
		   T * mean,
		   T * M2) {

  T delta = data - (*mean);
  (*mean) += delta / num;
  (*M2) += delta * (data - (*mean));
}

}

// ================================================================

#endif  // #ifndef RESULTS_ZENO_H_

