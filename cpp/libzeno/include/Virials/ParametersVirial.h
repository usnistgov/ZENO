// ================================================================
//
/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
//
// ================================================================

// ================================================================
// 
// Authors:
// Created:
//
// ================================================================

#ifndef PARAMETERS_VIRIALS_H_
#define PARAMETERS_VIRIALS_H_

// ================================================================

#include "Geometry/Vector3.h"
#include "Geometry/Sphere.h"

// ================================================================

namespace zeno {
  
/// Collects the parameters that are used when performing the Interior Sampling
/// computation.
///
/// For some parameters, tracks whether they 
/// have been manually set or are still at their default value.
///
class ParametersVirial
{
public:
  ParametersVirial();
  ~ParametersVirial();

  void setNumThreads(int numThreads);
  int getNumThreads() const;

  void setSeed(int seed);
  int getSeed() const;

  void setSteps(long long steps);
  long long getSteps() const;
  bool getStepsWasSet() const;

  void setMaxErrorVirialCoefficient(double maxErrorVirialCoefficient);
  double getMaxErrorVirialCoefficient() const;
  bool getMaxErrorVirialCoefficientWasSet() const;

  void setMaxRunTime(double maxRunTime);
  double getMaxRunTime() const;
  bool getMaxRunTimeWasSet() const;

  void setOrder(int order);
  int getOrder() const;

  void mpiBroadcast(int root);

private:
  void serializeMpiBroadcast(int root) const;
  void mpiBroadcastDeserialize(int root);
  
  // Command-line parameters
  
  int numThreads;

  int seed;

  long long steps;
  bool stepsWasSet;

  double maxErrorVirialCoefficient;
  bool maxErrorVirialCoefficientWasSet;

  double maxRunTime;
  bool maxRunTimeWasSet;

  int order;
};

}

// ================================================================

#endif  // #ifndef PARAMETERS_INTERIOR_SAMPLING_H_

