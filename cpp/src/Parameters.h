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
// Created: Tue Jan 05 16:39:23 2016 EDT
//
// ================================================================

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

// ================================================================

#include <string>
#include <fstream>

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

  void print(std::ofstream * csvOutputFile) const;

  std::string getInputFileName() const;

  std::string getCsvOutputFileName() const;
  bool getCsvOutputFileNameWasSet() const;
  
  int getMpiSize() const;
  void setMpiSize(int mpiSize);

  int getMpiRank() const;
  void setMpiRank(int mpiRank);

  int getNumThreads() const;

  int getSeed() const;

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

  double getMaxRunTime() const;
  bool getMaxRunTimeWasSet() const;

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

  std::string csvOutputFileName;
  bool csvOutputFileNameWasSet;

  int mpiSize;
  int mpiRank;

  int numThreads;

  int seed;

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

  double maxRunTime;
  bool maxRunTimeWasSet;

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

  template <typename T>
  void printScalarValue(std::string const & prettyName,
	   	        std::string const & csvName,
		        std::string const & units,
		        T property,
		        std::ofstream * csvOutputFile) const;

  template <typename T>
  void printVector3Value(std::string const & prettyName,
	   	         std::string const & csvName,
		         std::string const & units,
		         Vector3<T> const & property,
		         std::ofstream * csvOutputFile) const;
};

// ================================================================

template <typename T>
void
Parameters::printScalarValue(std::string const & prettyName,
	   	             std::string const & csvName,
		             std::string const & units,
		             T property,
		             std::ofstream * csvOutputFile) const {

  std::cout << prettyName << ": " << property << " " << units << std::endl;
  
  *csvOutputFile << csvName << ",units," << units
		 << std::endl
		 << csvName << ",value," << property
		 << std::endl;
}

template <typename T>
void
Parameters::printVector3Value(std::string const & prettyName,
	   	              std::string const & csvName,
		              std::string const & units,
		              Vector3<T> const & property,
		              std::ofstream * csvOutputFile) const {

  std::cout << prettyName << ": " << property << " " << units << std::endl;
  
  for (int i = 0; i < 3; ++i) {
    *csvOutputFile << csvName << "[" << i << "],"
		   << "units,"
		   << units
		   << std::endl
		   << csvName << "[" << i << "],"
		   << "value,"
		   << property.get(i)
		   << std::endl;
  }
}

// ================================================================

#endif  // #ifndef PARAMETERS_H_

