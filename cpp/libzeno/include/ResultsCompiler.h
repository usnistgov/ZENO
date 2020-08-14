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

#ifndef RESULTS_COMPILER_H_
#define RESULTS_COMPILER_H_

// ================================================================

#include <vector>
#include <array>
#include <cstdint>
#include <fstream>

#include "ResultsZeno.h"
#include "ResultsInterior.h"
#include "Results.h"
#include "Uncertain.h"
#include "ParametersResults.h"

#include "Geometry/Vector3.h"
#include "Geometry/Matrix3x3.h"
#include "Geometry/Sphere.h"

// ================================================================

namespace zeno {
  
/// Takes results from the Walk-on-Spheres and Interior Sampling computations
/// and derives physical quantities from them.
///
class ResultsCompiler {
public:
  ResultsCompiler(ParametersResults const & parameters);

  ~ResultsCompiler();

  void compile(ResultsZeno const * resultsZeno,
	       ResultsInterior const * resultsInterior,
	       bool computeForm,
	       Results * results);

private:
  using BigUInt = uint_fast64_t;

  static const unsigned int numFormFactors = 81;

  Result<Uncertain<double> > 
  computeCapacitance(Uncertain<double> const & t, 
		     double boundingSphereRadius) const;

  Result<Matrix3x3<Uncertain<double> > > 
  computePolarizability(Uncertain<double> const & t,
			Vector3<Uncertain<double> > const & u,
			Matrix3x3<Uncertain<double> > const & v,
			Matrix3x3<Uncertain<double> > const & w,
			double boundingSphereRadius) const;

  Result<Uncertain<double> > 
  computeMeanPolarizability(Matrix3x3<Uncertain<double> > const & 
			    polarizabilityTensor) const;

  Result<Uncertain<double> > 
  computePadeApproximant(Matrix3x3<Uncertain<double> > const & 
			 polarizabilityTensor) const;

  Uncertain<double> 
  computePadeApproximant(Uncertain<double> x1, 
			 Uncertain<double> x2) const;

  Result<Uncertain<double> >
  computeVolume(Uncertain<double> const & numInteriorHits, 
		double numInteriorSamples,
		double boundingSphereVolume) const;

  Result<Uncertain<double> >
  computeIntrinsicConductivity(Uncertain<double> const & meanPolarizability,
			       Uncertain<double> const & volume) const;

  Result<Uncertain<double> >
  computeCapacitanceOfASphere(Uncertain<double> const & volume) const;

  Result<Uncertain<double> >
  computeHydrodynamicRadius(Uncertain<double> const & capacitance) const;

  Result<Uncertain<double> > 
  computeViscometricRadius(Uncertain<double> const & meanPolarizability,
			   Uncertain<double> const & padeApproximant) const;

  Result<Uncertain<double> > 
  computeIntrinsicViscosity(Uncertain<double> const & padeApproximant,
			    Uncertain<double> const & intrinsicConductivity) 
    const;
			    
  Result<Uncertain<double> >
  computeIntrinsicViscosityConventional
  (Uncertain<double> const & padeApproximant,
   Uncertain<double> const & meanPolarizability,
   double mass) const;

  Result<Uncertain<double> >
  computeFrictionCoefficient(double solventViscosity,
			     Uncertain<double> const & hydrodynamicRadius) 
    const;

  Result<Uncertain<double> >
  computeDiffusionCoefficient(double temperature,
			      double solventViscosity,
			      Uncertain<double> const & hydrodynamicRadius) 
    const;

  Result<Uncertain<double> >
  computeSedimentationCoefficient(double mass,
				  double buoyancyFactor,
				  double solventViscosity,
				  Uncertain<double> const & hydrodynamicRadius)
    const;

  Result<Matrix3x3<Uncertain<double> > >
  computeGyrationTensor(Matrix3x3<Uncertain<double> > const & hitPointsSqrSum,
			Vector3<Uncertain<double> > const & hitPointsSum,
			Uncertain<double> const & numInteriorHits) const;

  std::array<double, numFormFactors>
  computeFormFactorQs(double boundingSphereRadius) const;

  std::array<double, numFormFactors>
  computeFormFactors(std::vector<Vector3<double> > const & interiorPoints,
		     std::array<double, 
		     numFormFactors> const & 
		     formFactorQs,
		     int numThreads) const;

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

  void printScalar(std::string const & prettyName,
		   std::string const & csvName,
		   std::string const & units,
		   Uncertain<double> const & property,
		   std::ofstream * csvOutputFile) const;

  void printVector3(std::string const & prettyName,
		    std::string const & csvName,
		    std::string const & units,
		    Vector3<Uncertain<double> > const & property,
		    std::ofstream * csvOutputFile) const;

  void printMatrix3x3(std::string const & prettyName,
		      std::string const & csvName,
		      std::string const & units,
		      Matrix3x3<Uncertain<double> > const & property,
		      std::ofstream * csvOutputFile) const;
  
  ParametersResults const * parameters;
};

}

// ================================================================

#endif  // #ifndef RESULTS_COMPILER_H_

