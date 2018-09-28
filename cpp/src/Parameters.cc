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
// Created: Tue Jan 05 16:39:56 2016 EDT
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
    csvOutputFileName(),
    csvOutputFileNameWasSet(false),
    mpiSize(1),
    mpiRank(0),
    numThreads(),
    seed(),
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
    maxRunTime(),
    maxRunTimeWasSet(false),
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
/// contains.  If a config file is specified on the command line, this will be
/// parsed as well.
///
/// Also computes default values for numThreads and seed if they were not set.
///
void 
Parameters::parseCommandLine(int argc, char **argv) {
  gengetopt_args_info args_info;
  cmdline_parser_params *parser_params;

  parser_params = cmdline_parser_params_create();

  //command line parse
  
  parser_params->initialize      = 1;
  parser_params->override        = 0;
  parser_params->check_required  = 0;
  parser_params->check_ambiguity = 0;

  int cmdlineParseResult =
    cmdline_parser_ext(argc, argv, &args_info, parser_params);
  
  if (cmdlineParseResult != 0) {
    exit(EXIT_FAILURE);
  }

  //config file parse

  if (args_info.config_file_given) {
    parser_params->initialize      = 0;
    parser_params->override        = 0;
    parser_params->check_required  = 0;
    parser_params->check_ambiguity = 0;
  
    int configParseResult =
      cmdline_parser_config_file(args_info.config_file_arg,
				 &args_info, parser_params);

    if (configParseResult != 0) {
      exit(EXIT_FAILURE);
    }
  }

  //check if required options are present

  int cmdlineCheckResult =
    cmdline_parser_required(&args_info, argv[0]);

  if (cmdlineCheckResult != 0) {
    exit(EXIT_FAILURE);
  }

  //copy arg values
  
  inputFileName = args_info.input_file_arg;

  if (args_info.csv_output_file_given) {
    csvOutputFileName = args_info.csv_output_file_arg;
    csvOutputFileNameWasSet = true;
  }
  
  minTotalNumWalks   = args_info.min_num_walks_arg;
  minTotalNumSamples = args_info.min_num_interior_samples_arg;

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

  if (args_info.max_run_time_given) {
    maxRunTime = args_info.max_run_time_arg;
    maxRunTimeWasSet = true;
  }

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
      std::cerr << "Error randomly setting seed for random number generator"
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

  free(parser_params);
}

/// Prints the parameters.  Most parameters are not printed if they have not
/// been set.
/// 
void
Parameters::print(std::ofstream * csvOutputFile) const {
  printScalarValue("Input file", "input_file", "",
		   inputFileName, csvOutputFile);

  printScalarValue("Number of nodes", "num_nodes", "",
		   mpiSize, csvOutputFile);

  printScalarValue("Number of threads", "num_threads", "",
		   numThreads, csvOutputFile);

  printScalarValue("Random number seed", "random_number_seed", "",
		   seed, csvOutputFile);

  if (totalNumWalksWasSet) {
    printScalarValue("Number of walks requested", "num_walks_requested", "",
		     totalNumWalks, csvOutputFile);
  }

  if (totalNumSamplesWasSet) {
    printScalarValue("Number of interior samples requested",
		     "num_interior_samples_requested", "",
		     totalNumSamples, csvOutputFile);
  }

  if (maxErrorCapacitanceWasSet) {
    printScalarValue("Max error in capacitance requested",
		     "max_error_capacitance_requested", "%",
	             maxErrorCapacitance, csvOutputFile);
  }

  if (maxErrorPolarizabilityWasSet) {
    printScalarValue("Max error in mean polarizability requested",
		     "max_error_mean_polarizability_requested", "%",
	             maxErrorPolarizability, csvOutputFile);
  }

  if (maxErrorVolumeWasSet) {
    printScalarValue("Max error in volume requested",
		     "max_error_volume_requested", "%",
	             maxErrorVolume, csvOutputFile);
  }

  if (maxRunTimeWasSet) {
    printScalarValue("Max run time", "max_run_time", "s",
		     maxRunTime, csvOutputFile);
  }

  if (skinThicknessWasSet) {
    printScalarValue("Skin thickness", "skin_thickness", "",
		     skinThickness, csvOutputFile);
  }

  if (launchCenterWasSet) {
    printVector3Value("Launch center", "launch_center", "",
		      launchCenter, csvOutputFile);
  }

  if (launchRadiusWasSet) {
    printScalarValue("Launch radius", "launch_radius", "",
		     launchRadius, csvOutputFile);
  }

  if (lengthScaleWasSet) {
    printScalarValue("Length scale", "length_scale",
		     Units::getName(lengthScaleUnit),
		     lengthScale, csvOutputFile);
  }

  if (temperatureWasSet) {
    printScalarValue("Temperature", "temperature",
		     Units::getName(temperatureUnit),
		     temperature, csvOutputFile);
  }

  if (massWasSet) {
    printScalarValue("Mass", "mass",
		     Units::getName(massUnit),
		     mass, csvOutputFile);
  }

  if (solventViscosityWasSet) {
    printScalarValue("Solvent viscosity", "solvent_viscosity",
		     Units::getName(solventViscosityUnit),
		     solventViscosity, csvOutputFile);
  }

  if (buoyancyFactorWasSet) {
    printScalarValue("Buoyancy factor", "buoyancy_factor", "",
		     buoyancyFactor, csvOutputFile);
  }
}

std::string 
Parameters::getInputFileName() const { 
  return inputFileName;
}

std::string
Parameters::getCsvOutputFileName() const {
  return csvOutputFileName;
}

bool
Parameters::getCsvOutputFileNameWasSet() const {
  return csvOutputFileNameWasSet;
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

double
Parameters::getMaxRunTime() const {
  return maxRunTime;
}

bool
Parameters::getMaxRunTimeWasSet() const {
  return maxRunTimeWasSet;
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
