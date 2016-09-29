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
// Time-stamp: <2016-08-29 14:42:11 dcj>
//
// ================================================================

#ifndef RESULTS_COMPILER_H_
#define RESULTS_COMPILER_H_

// ================================================================

#include <vector>
#include <array>
#include <cstdint>

#include "ResultsZeno.h"
#include "ResultsInterior.h"
#include "Uncertain.h"
#include "Parameters.h"

#include "Geometry/Vector3.h"
#include "Geometry/Matrix3x3.h"
#include "Geometry/Sphere.h"

// ================================================================

/// Takes results from the Walk-on-Spheres and Interior Sampling computations
/// and derives physical quantities from them.
///
class ResultsCompiler {
public:
  ResultsCompiler(Parameters const & parameters);

  ~ResultsCompiler();

  void compile(ResultsZeno const * resultsZeno,
	       ResultsInterior const * resultsInterior,
	       Sphere<double> const & boundingSphere,
	       bool computeForm);

  void print(bool printCounts) const;

  Uncertain<double> getCapacitance() const;

  Uncertain<double> getMeanPolarizability() const;

  Uncertain<double> getVolume() const;

private:
  using BigUInt = uint_fast64_t;

  static const unsigned int numFormFactors = 81;

  Uncertain<double> 
  computeCapacitance(Uncertain<double> const & t, 
		     double boundingSphereRadius) const;

  Matrix3x3<Uncertain<double> > 
  computePolarizability(Uncertain<double> const & t,
			Vector3<Uncertain<double> > const & u,
			Matrix3x3<Uncertain<double> > const & v,
			Matrix3x3<Uncertain<double> > const & w,
			double boundingSphereRadius) const;

  Uncertain<double> 
  computeMeanPolarizability(Matrix3x3<Uncertain<double> > const & 
			    polarizabilityTensor) const;

  Uncertain<double> 
  computePadeApproximant(Matrix3x3<Uncertain<double> > const & 
			 polarizabilityTensor) const;

  Uncertain<double> 
  computePadeApproximant(Uncertain<double> x1, 
			 Uncertain<double> x2) const;

  Uncertain<double>
  computeVolume(Uncertain<double> const & numInteriorHits, 
		double numInteriorSamples,
		double boundingSphereVolume) const;

  Uncertain<double>
  computeIntrinsicConductivity(Uncertain<double> const & meanPolarizability,
			       Uncertain<double> const & volume) const;

  Uncertain<double>
  computeCapacitanceOfASphere(Uncertain<double> const & volume) const;

  Uncertain<double>
  computeHydrodynamicRadius(Uncertain<double> const & capacitance) const;

  Uncertain<double> 
  computeViscometricRadius(Uncertain<double> const & meanPolarizability,
			   Uncertain<double> const & padeApproximant) const;

  Uncertain<double> 
  computeIntrinsicViscosity(Uncertain<double> const & padeApproximant,
			    Uncertain<double> const & intrinsicConductivity) 
    const;
			    
  Uncertain<double>
  computeIntrinsicViscosityConventional
  (Uncertain<double> const & padeApproximant,
   Uncertain<double> const & meanPolarizability,
   double mass) const;

  Uncertain<double>
  computeFrictionCoefficient(double solventViscosity,
			     Uncertain<double> const & hydrodynamicRadius) 
    const;

  Uncertain<double>
  computeDiffusionCoefficient(double temperature,
			      double solventViscosity,
			      Uncertain<double> const & hydrodynamicRadius) 
    const;

  Uncertain<double>
  computeSedimentationCoefficient(double mass,
				  double buoyancyFactor,
				  double solventViscosity,
				  Uncertain<double> const & hydrodynamicRadius)
    const;

  Matrix3x3<Uncertain<double> >
  computeGyrationTensor(Matrix3x3<Uncertain<double> > const & hitPointsSqrSum,
			Vector3<Uncertain<double> > const & hitPointsSum,
			Uncertain<double> const & numInteriorHits) const;

  std::array<double, numFormFactors>
  computeFormFactorQs(double boundingSphereRadius) const;

  std::array<double, numFormFactors>
  computeFormFactors(std::vector<Vector3<double> > const & interiorPoints,
		     std::array<double, 
		     numFormFactors> const & 
		     formFactorQs) const;

  void
  computeFormFactorsThread(int threadNum,
			   std::vector<Vector3<double> > const & 
			   interiorPoints,
			   std::array<double, 
			   numFormFactors> const & 
			   formFactorQs,
			   BigUInt startPairIndex,
			   BigUInt endPairIndex,
			   std::array<double, numFormFactors> * 
			   threadFormFactors) const;

  void indexToIJ(BigUInt index, BigUInt * i, BigUInt * j) const;

  Parameters const * parameters;

  double boundingSphereRadius;
  Vector3<double> boundingSphereCenter;

  Uncertain<double> t;
  Vector3<Uncertain<double> > u;
  Matrix3x3<Uncertain<double> > v;
  Matrix3x3<Uncertain<double> > w;

  Uncertain<double> numInteriorHits;

  Uncertain<double> capacitance;
  Matrix3x3<Uncertain<double> > polarizabilityTensor;
  Uncertain<double> meanPolarizability;
  Vector3<Uncertain<double> > polarizabilityEigenvalues;
  Uncertain<double> volume;
  Uncertain<double> intrinsicConductivity;
  Uncertain<double> capacitanceOfASphere;
  Uncertain<double> hydrodynamicRadius;
  Uncertain<double> q_eta;
  Uncertain<double> viscometricRadius;
  Uncertain<double> intrinsicViscosity;
  Uncertain<double> intrinsicViscosityConventional;
  Uncertain<double> frictionCoefficient;
  Uncertain<double> diffusionCoefficient;
  Uncertain<double> sedimentationCoefficient;
  Matrix3x3<Uncertain<double> > gyrationTensor;
  Vector3<Uncertain<double> > gyrationEigenvalues;
  std::array<double, numFormFactors> formFactorQs;
  std::array<double, numFormFactors> formFactors;

  bool intrinsicViscosityConventionalComputed;
  bool frictionCoefficientComputed;
  bool diffusionCoefficientComputed;
  bool sedimentationCoefficientComputed;

  bool resultsZenoCompiled;
  bool resultsInteriorCompiled;
  bool formResultsCompiled;
};

// ================================================================

#endif  // #ifndef RESULTS_COMPILER_H_

// ================================================================

// Local Variables:
// time-stamp-line-limit: 30
// mode: c++
// End:
