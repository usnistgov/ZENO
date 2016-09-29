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
// Date:    Tue Jan 05 16:39:23 2016 EDT
//
// Time-stamp: <2016-09-28 16:15:01 dcj>
//
// ================================================================

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

// ================================================================

#include <string>

#include "Units.h"
#include "Geometry/Vector3.h"

// ================================================================

/// Collects all the program parameters in one place.  Includes both
/// parameters read from the command line and from the bod file.
/// For some parameters, tracks whether they 
/// have been manually set or are still at their default value.
///
class Parameters
{
public:
  Parameters();
  ~Parameters();

  void parseCommandLine(int argc, char **argv);

  void print() const;

  std::string getInputFileName() const; 

  int getMpiSize() const;
  void setMpiSize(int mpiSize);

  int getMpiRank() const;
  void setMpiRank(int mpiRank);

  int getNumThreads() const;

  int getSeed() const;

  double getFracErrorBound() const;

  long long getTotalNumWalks() const;
  bool getTotalNumWalksWasSet() const;

  long long getTotalNumSamples() const;
  bool getTotalNumSamplesWasSet() const;

  double getMaxErrorCapacitance() const;
  bool getMaxErrorCapacitanceWasSet() const;

  double getMaxErrorPolarizability() const;
  bool getMaxErrorPolarizabilityWasSet() const;

  double getMaxErrorVolume() const;
  bool getMaxErrorVolumeWasSet() const;

  bool getComputeFormWasSet() const;

  long long getMinTotalNumWalks() const;

  long long getMinTotalNumSamples() const;

  std::string getSurfacePointsFileName() const;
  std::string getInteriorPointsFileName() const;

  bool getPrintCounts() const;
  bool getPrintBenchmarks() const;

  void setSkinThickness(double skinThickness);
  double getSkinThickness() const;
  bool getSkinThicknessWasSet() const;

  void setLaunchCenter(Vector3<double> launchCenter);
  Vector3<double> getLaunchCenter() const;
  bool getLaunchCenterWasSet() const;

  void setLaunchRadius(double launchRadius);
  double getLaunchRadius() const;
  bool getLaunchRadiusWasSet() const;

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

  void mpiSend() const;
  void mpiReceive();

private:
  std::string inputFileName; 

  int mpiSize;
  int mpiRank;

  int numThreads;

  int seed;

  double fracErrorBound;

  long long totalNumWalks;
  bool totalNumWalksWasSet;

  long long totalNumSamples;
  bool totalNumSamplesWasSet;

  double maxErrorCapacitance;
  bool maxErrorCapacitanceWasSet;

  double maxErrorPolarizability;
  bool maxErrorPolarizabilityWasSet;

  double maxErrorVolume;
  bool maxErrorVolumeWasSet;

  bool computeFormWasSet;

  long long minTotalNumWalks;

  long long minTotalNumSamples;

  std::string surfacePointsFileName;
  std::string interiorPointsFileName;

  bool printCounts;
  bool printBenchmarks;

  // .bod parameters

  double skinThickness;
  bool skinThicknessWasSet;

  Vector3<double> launchCenter;
  bool launchCenterWasSet;

  double launchRadius;
  bool launchRadiusWasSet;

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

// ================================================================

#endif  // #ifndef PARAMETERS_H_

// ================================================================

// Local Variables:
// time-stamp-line-limit: 30
// mode: c++
// End:
