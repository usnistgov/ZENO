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
// Author:  Derek Juba <derek.juba@nist.gov>
// Date:    Tue Jan 05 16:39:56 2016 EDT
// 
// Time-stamp: <2016-09-28 16:15:33 dcj>
// 
// ================================================================

#ifdef USE_MPI
#include <mpi.h>
#endif

#include <thread>
#include <cassert>
#include <iostream>
#include <fstream>

#include "cmdline.h"

#include "Parameters.h"

// ================================================================

Parameters::Parameters() 
  : inputFileName(), 
    mpiSize(1),
    mpiRank(0),
    numThreads(),
    seed(),
    fracErrorBound(),
    totalNumWalks(),
    totalNumWalksWasSet(false),
    totalNumSamples(),
    totalNumSamplesWasSet(false),
    maxErrorCapacitance(),
    maxErrorCapacitanceWasSet(false),
    maxErrorPolarizability(),
    maxErrorPolarizabilityWasSet(false),
    maxErrorVolume(),
    maxErrorVolumeWasSet(false),
    computeFormWasSet(false),
    minTotalNumWalks(),
    minTotalNumSamples(),
    surfacePointsFileName(),
    interiorPointsFileName(),
    printCounts(),
    printBenchmarks(),
    skinThickness(),
    skinThicknessWasSet(false),
    launchCenter(),
    launchCenterWasSet(false),
    launchRadius(),
    launchRadiusWasSet(false),
    lengthScale(),
    lengthScaleUnit(),
    lengthScaleWasSet(false),
    temperature(),
    temperatureUnit(),
    temperatureWasSet(false),
    mass(),
    massUnit(),
    massWasSet(false),
    solventViscosity(),
    solventViscosityUnit(),
    solventViscosityWasSet(false),
    buoyancyFactor(),
    buoyancyFactorWasSet(false) {

}

Parameters::~Parameters() {

}

/// Parses the command line and stores the values of the parameters it
/// contains.
///
/// Also sets default values for numThreads and seed.
///
void 
Parameters::parseCommandLine(int argc, char **argv) {
  struct gengetopt_args_info args_info;
  struct cmdline_parser_params *params;

  params = cmdline_parser_params_create();

  params->initialize = 1;
  params->check_required = 1;

  if (cmdline_parser_ext(argc, argv, &args_info, params) != 0) {
    cmdline_parser_free(&args_info);
    free(params);
    exit(1);
  }

  inputFileName      = args_info.input_file_arg;
  minTotalNumWalks   = args_info.min_num_walks_arg;
  minTotalNumSamples = args_info.min_num_interior_samples_arg;
  fracErrorBound     = args_info.frac_error_bound_arg;

  if (args_info.num_walks_given) {
    totalNumWalks = args_info.num_walks_arg;
    totalNumWalksWasSet = true;
  }

  if (args_info.num_interior_samples_given) {
    totalNumSamples = args_info.num_interior_samples_arg;
    totalNumSamplesWasSet = true;
  }

  if (args_info.max_rsd_capacitance_given) {
    maxErrorCapacitance = args_info.max_rsd_capacitance_arg;
    maxErrorCapacitanceWasSet = true;
  }

  if (args_info.max_rsd_polarizability_given) {
    maxErrorPolarizability = args_info.max_rsd_polarizability_arg;
    maxErrorPolarizabilityWasSet = true;
  }

  if (args_info.max_rsd_volume_given) {
    maxErrorVolume = args_info.max_rsd_volume_arg;
    maxErrorVolumeWasSet = true;
  }

  // if (args_info.compute_form_given) {
  //   computeFormWasSet = true;
  // }

  if (args_info.num_threads_given) {
    numThreads = args_info.num_threads_arg;
  }
  else {
    numThreads = std::thread::hardware_concurrency();

    if (numThreads == 0) {
      numThreads = 1;
    }
  }

  if (args_info.seed_given) {
    seed = args_info.seed_arg;
  }
  else {
    std::ifstream urandom("/dev/urandom");

    urandom.read((char *)&seed, sizeof(seed));

    if (!urandom.good()) {
      std::cout << "Error randomly setting seed for random number generator"
		<< std::endl;
    }

    urandom.close();

    seed = abs(seed);
  }

  if (args_info.surface_points_file_given) {
    surfacePointsFileName = args_info.surface_points_file_arg;
  }
  else {
    surfacePointsFileName = "";
  }

  if (args_info.interior_points_file_given) {
    interiorPointsFileName = args_info.interior_points_file_arg;
  }
  else {
    interiorPointsFileName = "";
  }

  printCounts     = args_info.print_counts_given;
  printBenchmarks = args_info.print_benchmarks_given;

  free(params);
}

/// Prints the parameters.  Most parameters are not printed if they have not
/// been set.
/// 
void
Parameters::print() const {
  std::cout << "Input file: " <<  inputFileName << std::endl 
	    << "Number of nodes: " << mpiSize << std::endl
	    << "Number of threads: " << numThreads << std::endl
	    << "Random number seed: " << seed << std::endl;

  if (totalNumWalksWasSet) {
    std::cout << "Number of walks: " << totalNumWalks << std::endl;
  }

  if (totalNumSamplesWasSet) {
    std::cout << "Number of interior samples: " << totalNumSamples << std::endl;
  }

  if (maxErrorCapacitanceWasSet) {
    std::cout << "Max error in capacitance: " << maxErrorCapacitance << " %" 
	      << std::endl;
  }

  if (maxErrorPolarizabilityWasSet) {
    std::cout << "Max error in mean polarizability: " << maxErrorPolarizability
	      << " %" << std::endl;
  }

  if (maxErrorVolumeWasSet) {
    std::cout << "Max error in volume: " << maxErrorVolume << " %" << std::endl;
  }

  if (skinThicknessWasSet) {
    std::cout << "Skin thickness: " << skinThickness << std::endl;
  }

  if (launchCenterWasSet) {
    std::cout << "Launch center: " << launchCenter << std::endl;
  }

  if (launchRadiusWasSet) {
    std::cout << "Launch radius: " << launchRadius << std::endl;
  }

  if (lengthScaleWasSet) {
    std::cout << "Length scale: " << lengthScale << " "
	      << Units::getName(lengthScaleUnit) << std::endl;
  }

  if (temperatureWasSet) {
    std::cout << "Temperature: " << temperature << " "
	      << Units::getName(temperatureUnit) << std::endl;
  }

  if (massWasSet) {
    std::cout << "Mass: " << mass << " "
	      << Units::getName(massUnit) << std::endl;
  }

  if (solventViscosityWasSet) {
    std::cout << "Solvent viscosity: " << solventViscosity << " "
	      << Units::getName(solventViscosityUnit) << std::endl;
  }

  if (buoyancyFactorWasSet) {
    std::cout << "Buoyancy factor: " << buoyancyFactor << std::endl;
  }
}

std::string 
Parameters::getInputFileName() const { 
  return inputFileName;
}

int
Parameters::getMpiSize() const {
  return mpiSize;
}

void
Parameters::setMpiSize(int mpiSize) {
  this->mpiSize = mpiSize;
}

int
Parameters::getMpiRank() const {
  return mpiRank;
}

void
Parameters::setMpiRank(int mpiRank) {
  this->mpiRank = mpiRank;
}

int 
Parameters::getNumThreads() const {
  return numThreads;
}

int 
Parameters::getSeed() const {
  return seed;
}

double 
Parameters::getFracErrorBound() const {
  return fracErrorBound;
}

long long 
Parameters::getTotalNumWalks() const {
  return totalNumWalks;
}

bool
Parameters::getTotalNumWalksWasSet() const {
  return totalNumWalksWasSet;
}

long long 
Parameters::getTotalNumSamples() const {
  return totalNumSamples;
}

bool 
Parameters::getTotalNumSamplesWasSet() const {
  return totalNumSamplesWasSet;
}

double
Parameters::getMaxErrorCapacitance() const {
  return maxErrorCapacitance;
}

bool
Parameters::getMaxErrorCapacitanceWasSet() const {
  return maxErrorCapacitanceWasSet;
}

double
Parameters::getMaxErrorPolarizability() const {
  return maxErrorPolarizability;
}

bool
Parameters::getMaxErrorPolarizabilityWasSet() const {
  return maxErrorPolarizabilityWasSet;
}

bool
Parameters::getComputeFormWasSet() const {
  return computeFormWasSet;
}

double
Parameters::getMaxErrorVolume() const {
  return maxErrorVolume;
}

bool
Parameters::getMaxErrorVolumeWasSet() const {
  return maxErrorVolumeWasSet;
}

long long 
Parameters::getMinTotalNumWalks() const {
  return minTotalNumWalks;
}

long long 
Parameters::getMinTotalNumSamples() const {
  return minTotalNumSamples;
}

std::string 
Parameters::getSurfacePointsFileName() const {
  return surfacePointsFileName;
}

std::string 
Parameters::getInteriorPointsFileName() const {
  return interiorPointsFileName;
}

bool 
Parameters::getPrintCounts() const {
  return printCounts;
}

bool 
Parameters::getPrintBenchmarks() const {
  return printBenchmarks;
}

void 
Parameters::setSkinThickness(double skinThickness) {
  this->skinThickness = skinThickness;

  skinThicknessWasSet = true;
}

double
Parameters::getSkinThickness() const {
  return skinThickness;
}

bool
Parameters::getSkinThicknessWasSet() const {
  return skinThicknessWasSet;
}

void 
Parameters::setLaunchCenter(Vector3<double> launchCenter) {
  this->launchCenter = launchCenter;

  launchCenterWasSet = true;
}

Vector3<double>
Parameters::getLaunchCenter() const {
  return launchCenter;
}

bool 
Parameters::getLaunchCenterWasSet() const {
  return launchCenterWasSet;
}

void 
Parameters::setLaunchRadius(double launchRadius) {
  this->launchRadius = launchRadius;

  launchRadiusWasSet = true;
}

double
Parameters::getLaunchRadius() const {
  return launchRadius;
}

bool 
Parameters::getLaunchRadiusWasSet() const {
  return launchRadiusWasSet;
}

void
Parameters::setLengthScale(double number, Units::Length unit) {
  lengthScale     = number;
  lengthScaleUnit = unit;

  lengthScaleWasSet = true;
}

double
Parameters::getLengthScaleNumber() const {
  return lengthScale;
}

Units::Length
Parameters::getLengthScaleUnit() const {
  return lengthScaleUnit;
}

bool
Parameters::getLengthScaleWasSet() const {
  return lengthScaleWasSet;
}

void 
Parameters::setTemperature(double number, Units::Temperature unit) {
  temperature     = number;
  temperatureUnit = unit;

  temperatureWasSet = true;
}

double
Parameters::getTemperatureNumber() const {
  return temperature;
}

Units::Temperature
Parameters::getTemperatureUnit() const {
  return temperatureUnit;
}

bool
Parameters::getTemperatureWasSet() const {
  return temperatureWasSet;
}

void
Parameters::setMass(double number, Units::Mass unit) {
  mass     = number;
  massUnit = unit;

  massWasSet = true;
}

double 
Parameters::getMassNumber() const {
  return mass;
}

Units::Mass
Parameters::getMassUnit() const {
  return massUnit;
}

bool
Parameters::getMassWasSet() const {
  return massWasSet;
}

void 
Parameters::setSolventViscosity(double number, Units::Viscosity unit) {
  solventViscosity     = number;
  solventViscosityUnit = unit;

  solventViscosityWasSet = true;
}

double
Parameters::getSolventViscosityNumber() const {
  return solventViscosity;
}

Units::Viscosity
Parameters::getSolventViscosityUnit() const {
  return solventViscosityUnit;
}

bool
Parameters::getSolventViscosityWasSet() const {
  return solventViscosityWasSet;
}

void 
Parameters::setBuoyancyFactor(double buoyancyFactor) {
  this->buoyancyFactor = buoyancyFactor;

  buoyancyFactorWasSet = true;
}

double
Parameters::getBuoyancyFactor() const {
  return buoyancyFactor;
}

bool
Parameters::getBuoyancyFactorWasSet() const {
  return buoyancyFactorWasSet;
}

/// Broadcasts the parameters that can be set in the bod file over MPI.
///
void 
Parameters::mpiSend() const {
#ifdef USE_MPI
  double parametersArray[22];

  parametersArray[0]  = getSkinThickness();
  parametersArray[1]  = (double)getSkinThicknessWasSet();
  parametersArray[2]  = getLaunchCenter().getX();
  parametersArray[3]  = getLaunchCenter().getY();
  parametersArray[4]  = getLaunchCenter().getZ();
  parametersArray[5]  = (double)getLaunchCenterWasSet();
  parametersArray[6]  = getLaunchRadius();
  parametersArray[7]  = (double)getLaunchRadiusWasSet();
  parametersArray[8]  = getLengthScaleNumber();
  parametersArray[9]  = (double)getLengthScaleUnit();
  parametersArray[10] = (double)getLengthScaleWasSet();
  parametersArray[11] = getTemperatureNumber();
  parametersArray[12] = (double)getTemperatureUnit();
  parametersArray[13] = (double)getTemperatureWasSet();
  parametersArray[14] = getMassNumber();
  parametersArray[15] = (double)getMassUnit();
  parametersArray[16] = (double)getMassWasSet();
  parametersArray[17] = getSolventViscosityNumber();
  parametersArray[18] = (double)getSolventViscosityUnit();
  parametersArray[19] = (double)getSolventViscosityWasSet();
  parametersArray[20] = getBuoyancyFactor();
  parametersArray[21] = (double)getBuoyancyFactorWasSet();

  MPI_Bcast(parametersArray, 22, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
}

/// Receives the parameters that can be set in the bod file over MPI.
///
void 
Parameters::mpiReceive() {
#ifdef USE_MPI
  double parametersArray[22];

  MPI_Bcast(parametersArray, 22, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if ((bool)parametersArray[1]) {
    setSkinThickness(parametersArray[0]);
  }

  if ((bool)parametersArray[5]) {
    setLaunchCenter(Vector3<double>(parametersArray[2],
				    parametersArray[3],
				    parametersArray[4]));
  }

  if ((bool)parametersArray[7]) {
    setLaunchRadius(parametersArray[6]);
  }

  if ((bool)parametersArray[10]) {
    setLengthScale(parametersArray[8], 
		   (Units::Length)parametersArray[9]);
  }

  if ((bool)parametersArray[13]) {
    setTemperature(parametersArray[11],
		   (Units::Temperature)parametersArray[12]);
  }

  if ((bool)parametersArray[16]) {
    setMass(parametersArray[14],
	    (Units::Mass)parametersArray[15]);
  }

  if ((bool)parametersArray[19]) {
    setSolventViscosity(parametersArray[17],
			(Units::Viscosity)parametersArray[18]);
  }

  if ((bool)parametersArray[21]) {
    setBuoyancyFactor(parametersArray[20]);
  }
#endif
}

// ================================================================

// Local Variables:
// time-stamp-line-limit: 30
// End:
