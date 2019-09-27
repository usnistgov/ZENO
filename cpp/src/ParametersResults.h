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

#ifndef PARAMETERS_RESULTS_H_
#define PARAMETERS_RESULTS_H_

// ================================================================

#include "Units.h"

// ================================================================

namespace zeno {
  
/// Collects the parameters that are used when computing results.
///
/// For some parameters, tracks whether they 
/// have been manually set or are still at their default value.
///
class ParametersResults
{
public:
  ParametersResults();
  ~ParametersResults();

  void setComputeForm(bool computeForm);
  bool getComputeForm() const;
  
  void setLengthScale(double number, Units::Length unit);
  double getLengthScaleNumber() const;
  Units::Length getLengthScaleUnit() const;
  bool getLengthScaleWasSet() const;

  void setTemperature(double number, Units::Temperature unit);
  double getTemperatureNumber() const;
  Units::Temperature getTemperatureUnit() const;
  bool getTemperatureWasSet() const;

  void setMass(double number, Units::Mass unit);
  double getMassNumber() const;
  Units::Mass getMassUnit() const;
  bool getMassWasSet() const;

  void setSolventViscosity(double number, Units::Viscosity unit);
  double getSolventViscosityNumber() const;
  Units::Viscosity getSolventViscosityUnit() const;
  bool getSolventViscosityWasSet() const;

  void setBuoyancyFactor(double buoyancyFactor);
  double getBuoyancyFactor() const;
  bool getBuoyancyFactorWasSet() const;

  void mpiBroadcast(int root);

private:
  void serializeMpiBroadcast(int root) const;
  void mpiBroadcastDeserialize(int root);
  
  // Command-line parameters

  bool computeForm;

  // .bod parameters
  
  double lengthScale;
  Units::Length lengthScaleUnit;
  bool lengthScaleWasSet;

  double temperature;
  Units::Temperature temperatureUnit;
  bool temperatureWasSet;

  double mass;
  Units::Mass massUnit;
  bool massWasSet;

  double solventViscosity;
  Units::Viscosity solventViscosityUnit;
  bool solventViscosityWasSet;

  double buoyancyFactor;
  bool buoyancyFactorWasSet;
};

}

// ================================================================

#endif  // #ifndef PARAMETERS_RESULTS_H_

