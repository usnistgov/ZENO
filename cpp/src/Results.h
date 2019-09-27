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
// Created: 2019-05-10
//
// ================================================================

#ifndef RESULTS_H
#define RESULTS_H

#include <array>

#include "Uncertain.h"

#include "Geometry/Vector3.h"
#include "Geometry/Matrix3x3.h"

#include "Result.h"

namespace zeno {

/// A struct collecting the results that can be computed by Zeno
///
class Results {
 public:
  Results();
  ~Results();

  static const unsigned int numFormFactors = 81;

  Result<long long> numWalks;
  
  Result<Uncertain<double> > t;
  Result<Vector3<Uncertain<double> > > u;
  Result<Matrix3x3<Uncertain<double> > > v;
  Result<Matrix3x3<Uncertain<double> > > w;

  Result<long long> numInteriorSamples;
  
  Result<Uncertain<double> > numInteriorHits;
  
  Result<Uncertain<double> > capacitance;
  Result<Matrix3x3<Uncertain<double> > > polarizabilityTensor;
  Result<Uncertain<double> > meanPolarizability;
  Result<Vector3<Uncertain<double> > > polarizabilityEigenvalues;
  Result<Uncertain<double> > volume;
  Result<Uncertain<double> > intrinsicConductivity;
  Result<Uncertain<double> > capacitanceOfASphere;
  Result<Uncertain<double> > hydrodynamicRadius;
  Result<Uncertain<double> > q_eta;
  Result<Uncertain<double> > viscometricRadius;
  Result<Uncertain<double> > intrinsicViscosity;
  Result<Uncertain<double> > intrinsicViscosityConventional;
  Result<Uncertain<double> > frictionCoefficient;
  Result<Uncertain<double> > diffusionCoefficient;
  Result<Uncertain<double> > sedimentationCoefficient;
  Result<Matrix3x3<Uncertain<double> > > gyrationTensor;
  Result<Vector3<Uncertain<double> > > gyrationEigenvalues;

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

}

#endif
