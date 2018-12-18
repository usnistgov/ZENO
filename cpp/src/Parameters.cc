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
#include <string>
#include <vector>

#include "Parameters.h"

#include <boost/program_options.hpp>

using std::string;
using std::vector;

namespace po = boost::program_options;

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
	
	// argument descriptions
	string configFileDesc = "Config file name";
	string inputFileDesc = "Input file name";
	string csvOutputFileDesc = "Write output in CSV format to the specified file in addition to displaying the regular output";
	string totalNumWalksDesc = "Number of walk-on-spheres walks to perform.  Fewer walks may be performed if another stopping condition is reached first";
	string totalNumSamplesDesc = "Number of interior samples to take.  Fewer samples may be taken if another stopping condition is reached first";
	string maxErrorCapacitanceDesc = "Perform walk-on-spheres walks until the relative standard deviation of the capacitance drops below this value.  This value may not be reached if another stopping condition is reached first.  Relative standard deviation is defined as (Standard_Deviation/Mean)*100%";
	string maxErrorPolarizabilityDesc = "Perform walk-on-spheres walks until the relative standard deviation of the mean electric polarizability drops below this value. This value may not be reached if another stopping condition is reached first. Relative standard deviation is defined as (Standard_Deviation/Mean)*100%";
	string maxErrorVolumeDesc = "Take interior samples until the relative standard deviation of volume drops below this value.  This value may not be reached if another stopping condition is reached first.  Relative standard deviation is defined as (Standard_Deviation/Mean)*100%";
	string minTotalNumWalksDesc = "Minimum number of walk-on-spheres walks to perform when using max-rsd stopping conditions";
	string minTotalNumSamplesDesc = "Minimum number of interior samples to take when using max-rsd stopping conditions";
	string maxRunTimeDesc = "Max time (in seconds) that the program is allowed to run.  If this time is reached, the computation will be stopped and the results computed so far will be displayed";
	string numThreadsDesc = "Number of threads to use  (default=Number of logical cores)";
	string seedDesc = "Seed for the random number generator  (default=Randomly set)";
	string surfacePointsFileDesc = "Name of file for writing the surface points from Walk-on-Spheres";
	string interiorPointsFileDesc = "Name of file for writing the interior sample points";
	string printCountsDesc = "Print statistics related to counts of hit points";
	string printBenchmarksDesc = "Print detailed RAM and timing information";
	
	po::options_description opts("zeno v5.1");
	opts.add_options()
		("help,h", "Print help and exit")
		("version,V", "Print version and exit")
		("config-file,c", po::value<string>(), configFileDesc.c_str())
		("input-file,i", po::value<string>(&inputFileName), inputFileDesc.c_str())
		("csv-output-file", po::value<string>(&csvOutputFileName), csvOutputFileDesc.c_str())
		("num-walks", po::value<long long>(&totalNumWalks), totalNumWalksDesc.c_str())
		("num-interior-samples", po::value<long long>(&totalNumSamples), totalNumSamplesDesc.c_str())
		("max-rsd-capacitance", po::value<double>(&maxErrorCapacitance), maxErrorCapacitanceDesc.c_str())
		("max-rsd-polarizability", po::value<double>(&maxErrorPolarizability), maxErrorPolarizabilityDesc.c_str())
		("max-rsd-volume", po::value<double>(&maxErrorVolume), maxErrorVolumeDesc.c_str())
		("min-num-walks", po::value<long long>(&minTotalNumWalks)->default_value(1000), minTotalNumWalksDesc.c_str())
		("min-num-interior-samples", po::value<long long>(&minTotalNumSamples)->default_value(1000), minTotalNumSamplesDesc.c_str())
		("max-run-time", po::value<double>(&maxRunTime), maxRunTimeDesc.c_str())
		("num-threads", po::value<int>(&numThreads), numThreadsDesc.c_str())
		("seed", po::value<int>(&seed), seedDesc.c_str())
		("surface-points-file", po::value<string>(&surfacePointsFileName), surfacePointsFileDesc.c_str())
		("interior-points-file", po::value<string>(&interiorPointsFileName), interiorPointsFileDesc.c_str())
		("print-counts", printCountsDesc.c_str())
		("print-benchmarks", printBenchmarksDesc.c_str());

	// argument map
	po::variables_map arg_map;
	// attempt to parse/store and exit on fail
	try {
		po::store(po::command_line_parser(argc, argv).options(opts).allow_unregistered().run(), arg_map);
		po::notify(arg_map);
	}
	catch(const std::exception& e) {
		std::cerr << e.what() << "\n";
		exit(EXIT_FAILURE);
	}
	
	// check for config file and parse if true
	if (arg_map.count("config-file")) {
		try {
			string configFileName = arg_map["config-file"].as<string>();
			po::store(po::parse_config_file<char>(configFileName.c_str(), opts), arg_map);
			po::notify(arg_map);
		}
		catch(const std::exception& e) {
			std::cerr << e.what() << "\n";
			exit(EXIT_FAILURE);
		}
	}
	
	// csv output file name
	if (arg_map.count("csv-output-file"))
		csvOutputFileNameWasSet = true;
	// number of walks
	if (arg_map.count("num-walks"))
		totalNumWalksWasSet = true;
	// number of interior samples
	if (arg_map.count("num-interior-samples"))
		totalNumSamplesWasSet = true;
	// max error capacitance
	if (arg_map.count("max-rsd-capacitance"))
		maxErrorCapacitanceWasSet = true;
	// max error polarizability
	if (arg_map.count("max-rsd-polarizability"))
		maxErrorPolarizabilityWasSet = true;
	// max error volume
	if (arg_map.count("max-rsd-volume"))
		maxErrorVolumeWasSet = true;
	// max run time
	if (arg_map.count("max-run-time"))
		maxRunTimeWasSet = true;
	// number of threads
	if (!arg_map.count("num-threads")) {
		numThreads = std::thread::hardware_concurrency();

		if (numThreads == 0)
			numThreads = 1;
	}
	// seed
	if (!arg_map.count("seed")) {
		std::ifstream urandom("/dev/urandom");

		urandom.read((char *)&seed, sizeof(seed));

		if (!urandom.good())
		  std::cerr << "Error randomly setting seed for random number generator" << std::endl;

		urandom.close();
		seed = abs(seed);
	}
	// surface points file
	if (!arg_map.count("surface-points-file"))
		surfacePointsFileName = "";
	// interior points file
	if (!arg_map.count("interior-points-file"))
		interiorPointsFileName = "";
	// print counts
	printCounts = false;
	if (arg_map.count("print-counts"))
		printCounts = true;
	// print benchmarks
	printBenchmarks = false;
	if (arg_map.count("print-benchmarks"))
		printBenchmarks = true;
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
