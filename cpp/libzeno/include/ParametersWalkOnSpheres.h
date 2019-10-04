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

#ifndef PARAMETERS_WALK_ON_SPHERES_H_
#define PARAMETERS_WALK_ON_SPHERES_H_

// ================================================================

#include "Geometry/Vector3.h"
#include "Geometry/Sphere.h"

// ================================================================

namespace zeno {
  
/// Collects the parameters that are used when performing the Walk-on-Spheres
/// computation.
///
/// For some parameters, tracks whether they 
/// have been manually set or are still at their default value.
///
class ParametersWalkOnSpheres
{
public:
  ParametersWalkOnSpheres();
  ~ParametersWalkOnSpheres();

  void setNumThreads(int numThreads);
  int getNumThreads() const;

  void setSeed(int seed);
  int getSeed() const;

  void setTotalNumWalks(long long totalNumWalks);
  long long getTotalNumWalks() const;
  bool getTotalNumWalksWasSet() const;

  void setMaxErrorCapacitance(double maxErrorCapacitance);
  double getMaxErrorCapacitance() const;
  bool getMaxErrorCapacitanceWasSet() const;

  void setMaxErrorPolarizability(double maxErrorPolarizability);
  double getMaxErrorPolarizability() const;
  bool getMaxErrorPolarizabilityWasSet() const;

  void setMaxRunTime(double maxRunTime);
  double getMaxRunTime() const;
  bool getMaxRunTimeWasSet() const;

  void setMinTotalNumWalks(long long minTotalNumWalks);
  long long getMinTotalNumWalks() const;

  void setSaveSurfacePoints(bool saveSurfacePoints);
  bool getSaveSurfacePoints() const;
  
  void setSkinThickness(double skinThickness);
  double getSkinThickness() const;
  bool getSkinThicknessWasSet() const;

  void setLaunchCenter(Vector3<double> launchCenter);
  Vector3<double> getLaunchCenter() const;
  bool getLaunchCenterWasSet() const;

  void setLaunchRadius(double launchRadius);
  double getLaunchRadius() const;
  bool getLaunchRadiusWasSet() const;

  Sphere<double> getLaunchSphere() const;

  void mpiBroadcast(int root);

private:
  void serializeMpiBroadcast(int root) const;
  void mpiBroadcastDeserialize(int root);
  
  // Command-line parameters
  
  int numThreads;

  int seed;

  long long totalNumWalks;
  bool totalNumWalksWasSet;

  double maxErrorCapacitance;
  bool maxErrorCapacitanceWasSet;

  double maxErrorPolarizability;
  bool maxErrorPolarizabilityWasSet;

  double maxRunTime;
  bool maxRunTimeWasSet;

  long long minTotalNumWalks;

  bool saveSurfacePoints;

  // .bod parameters

  double skinThickness;
  bool skinThicknessWasSet;

  Vector3<double> launchCenter;
  bool launchCenterWasSet;

  double launchRadius;
  bool launchRadiusWasSet;

};

}

// ================================================================

#endif  // #ifndef PARAMETERS_WALK_ON_SPHERES_H_
