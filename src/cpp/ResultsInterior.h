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
// Time-stamp: <2016-08-29 14:59:25 dcj>
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

// ================================================================

// Local Variables:
// time-stamp-line-limit: 30
// mode: c++
// End:
