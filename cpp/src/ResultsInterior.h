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

#ifndef RESULTS_INTERIOR_H_
#define RESULTS_INTERIOR_H_

// ================================================================

#include <vector>

#include "Geometry/Vector3.h"

#include "Uncertain.h"

// ================================================================

/// Collects results from the Interior Sampling computation.
///
class ResultsInterior {
public:
  ResultsInterior(int numThreads,
		  bool saveHitPoints);

  ~ResultsInterior();

  void recordHit(int threadNum,
		 Vector3<double> const & point);

  void recordMiss(int threadNum);

  void reduce();

  void gatherHitPoints();

  Uncertain<double> getNumHits() const;

  Matrix3x3<Uncertain<double> > getHitPointsSqrSum() const;

  Vector3<Uncertain<double> > getHitPointsSum() const;

  double getNumSamples() const;

  bool getSaveHitPoints() const;

  std::vector<Vector3<double> > const * getPoints() const;

private:
  void updateVariance(int threadNum,
		      double hitMissData,
		      Matrix3x3<double> const & hitPointsSqrData,
		      Vector3<double> const & hitPointsData);

  template <class T>
  void updateItemVariance(T const & data,
			  double num,
			  T * mean,
			  T * M2);

  template <class T>
  void reduceItem(T const & mean,
		  T const & M2,
		  double num,
		  T * sumReduced,
		  T * sumVarianceReduced);

  int const numThreads;

  bool saveHitPoints;

  double * numSamples;

  double * hitMissMean;
  double * hitMissM2;

  Matrix3x3<double> * hitPointsSqrMean;
  Matrix3x3<double> * hitPointsSqrM2;

  Vector3<double> * hitPointsMean;
  Vector3<double> * hitPointsM2;

  double numSamplesReduced;

  double numHitsReduced;
  double numHitsVarianceReduced;

  Matrix3x3<double> hitPointsSqrSumReduced;
  Matrix3x3<double> hitPointsSqrSumVarianceReduced;

  Vector3<double> hitPointsSumReduced;
  Vector3<double> hitPointsSumVarianceReduced;

  std::vector<Vector3<double> > * points;

  std::vector<Vector3<double> > gatheredPoints;

  bool reduced;
  bool hitPointsGathered;
};

template <class T>
void
ResultsInterior::
updateItemVariance(T const & data,
		   double num,
		   T * mean,
		   T * M2) {

  T delta = data - (*mean);
  (*mean) += delta / num;
  (*M2) += delta * (data - (*mean));
}

template <class T>
void
ResultsInterior::
reduceItem(T const & mean,
	   T const & M2,
	   double num,
	   T * sumReduced,
	   T * sumVarianceReduced) {

    T sum = mean * num;

    (*sumReduced) += sum;

    T variance = M2 / (num - 1);
    T meanVariance = variance / num;
    T sumVariance = meanVariance * pow(num, 2);

    (*sumVarianceReduced) += sumVariance;
}

// ================================================================

#endif  // #ifndef RESULTS_INTERIOR_H_

