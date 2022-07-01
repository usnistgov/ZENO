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
// Created: Mon Feb 24 11:25:18 2014 EDT
//
// ================================================================

#ifdef USE_MPI
#include <mpi.h>
#endif

#include <thread>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <limits>
#include <algorithm>
#include <array>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <ctime>

#include "cmdline.h"

#include "ParametersWalkOnSpheres.h"
#include "ParametersInteriorSampling.h"
#include "ParametersResults.h"
#include "ParametersLocal.h"

#include "Units.h"
#include "Uncertain.h"

#include "BodParser/BodParser.h"
#include "MapParser/MapParser.h"
#include "XyzParser/XyzParser.h"

#include "Results.h"

#include "Geometry/Vector3.h"
#include "Geometry/MixedModel.h"

#include "Timer.h"

#include "Zeno.h"

using namespace zeno;

// ================================================================

// [0] header
// [1] type
// [2] data
using CsvItems = std::array<std::vector<std::string>, 3>;

// ================================================================

void 
parseCommandLine(int argc, char **argv,
		 ParametersWalkOnSpheres * parametersWalkOnSpheres,
		 ParametersInteriorSampling * parametersInteriorSampling,
		 ParametersVirial * parametersVirial,
		 ParametersResults * parametersResults,
		 ParametersLocal * parametersLocal);

void
parseBodFile(ParametersLocal * parametersLocal,
	     ParametersWalkOnSpheres * parametersWalkOnSpheres,
	     ParametersInteriorSampling * parametersInteriorSampling,
	     ParametersResults * parametersResults,
	     MixedModel<double> * model);

void
parseMapFile(ParametersLocal const & parametersLocal,
	     std::unordered_map<std::string, double> * atomIdToRadius);

void
parseXyzFile(ParametersLocal const & parametersLocal,
	     std::unordered_map<std::string, double> const & atomIdToRadius,
	     std::list<zeno::MixedModel<double> > * snapshots,
	     std::list<std::string> * comments);

void
broadcastSnapshots(std::list<zeno::MixedModel<double> > * snapshots);
  
int
runZeno(ParametersLocal const & parametersLocal,
	ParametersWalkOnSpheres * parametersWalkOnSpheres,
	ParametersInteriorSampling * parametersInteriorSampling,
	ParametersVirial * parametersVirial,
	ParametersResults * parametersResults,
	double readTime,
	double broadcastTime,
	MixedModel<double> * model,
	CsvItems * csvItems);

void
printOutput(Results const & results,
	    ParametersWalkOnSpheres const & parametersWalkOnSpheres,
	    ParametersInteriorSampling const & parametersInteriorSampling,
	    ParametersVirial const & parametersVirial,
	    ParametersResults const & parametersResults,
	    ParametersLocal const & parametersLocal,
	    double initializeTime,
	    double readTime,
	    double broadcastTime,
	    double preprocessTime,
	    double walkTime,
	    double reduceTime,
	    double sampleTime,
	    double volumeReduceTime,
	    double virialTime,
	    double virialReduceTime,
	    double totalZenoTime,
	    CsvItems * csvItems);

void
printParameters(ParametersWalkOnSpheres const & parametersWalkOnSpheres,
	        ParametersInteriorSampling const & parametersInteriorSampling,
	        ParametersVirial const & parametersVirial,
	        ParametersResults const & parametersResults,
		ParametersLocal const & parametersLocal,
		CsvItems * csvItems);

void
printResults(Results const & results,
	     CsvItems * csvItems);

template <typename T>
void
printExactScalar(std::string const & prettyName,
		 std::string const & csvName,
		 std::string const & units,
		 T property,
		 CsvItems * csvItems);

template <typename T>
void
printExactVector3(std::string const & prettyName,
		  std::string const & csvName,
		  std::string const & units,
		  Vector3<T> const & property,
		  CsvItems * csvItems);

void
printScalar(Result<Uncertain<double> > const & result,
	    CsvItems * csvItems);

void
printVector3(Result<Vector3<Uncertain<double> > > const & result,
	     CsvItems * csvItems);

void
printMatrix3x3(Result<Matrix3x3<Uncertain<double> > > const & result,
	       CsvItems * csvItems);

void
savePointFiles(Zeno * zeno,
	       ParametersLocal const & parametersLocal);

void
printTime(std::string const & label);

void
printRAM(std::string const & label,
	 CsvItems * csvItems);

void
writeCsvOutputCols(CsvItems const & globalCsvItems,
		   std::list<CsvItems> const & perRunCsvItemsList,
		   std::ofstream * csvOutputFile);

void
writeCsvOutputRows(CsvItems const & globalCsvItems,
		   std::list<CsvItems> const & perRunCsvItemsList,
		   std::ofstream * csvOutputFile);

void
writePoints(std::string const & fileName, 
	    std::vector<Vector3<double> > const * points,
	    std::vector<Vector3<char> > const * charges);

template <typename T>
int
getRandomNumberFromOS(T * randomNumber);

template <typename T>
std::string
to_string_scientific(T const & value);

// ================================================================

int main(int argc, char **argv) {

  int mpiSize = 1, mpiRank = 0;

#ifdef USE_MPI
  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
#endif

  ParametersWalkOnSpheres parametersWalkOnSpheres;
  ParametersInteriorSampling parametersInteriorSampling;
  ParametersVirial parametersVirial;
  ParametersResults parametersResults;
  ParametersLocal parametersLocal;

  parametersLocal.setMpiSize(mpiSize);
  parametersLocal.setMpiRank(mpiRank);

  parseCommandLine(argc, argv,
		   &parametersWalkOnSpheres,
		   &parametersInteriorSampling,
		   &parametersVirial,
		   &parametersResults,
		   &parametersLocal);

  if (mpiRank == 0) {
    printTime("Start time: ");
  }

  std::ofstream csvOutputFile;

  if (parametersLocal.getCsvOutputFileNameWasSet() &&
      parametersLocal.getMpiRank() == 0) {
    
    csvOutputFile.open(parametersLocal.getCsvOutputFileName(),
		       std::ios::out | std::ios::trunc);

    if (!csvOutputFile.is_open()) {
      std::cerr << "*** Warning ***" << std::endl
		<< "Could not open CSV output file "
		<< parametersLocal.getCsvOutputFileName() << std::endl
		<< std::endl;
    }
  }

  CsvItems globalCsvItems;
  std::list<CsvItems> perRunCsvItemsList;

  if (parametersLocal.getMpiRank() == 0) {
    printExactScalar("Version",
		     "version",
		     "",
		     CMDLINE_PARSER_VERSION,
		     &globalCsvItems);
  }

  if (parametersLocal.getPrintBenchmarks() && 
      parametersLocal.getMpiRank() == 0) {

    printRAM("initialization",
	     &globalCsvItems);
  }

  MixedModel<double> model;

  Timer readTimer;

  if (parametersLocal.getMpiRank() == 0) {
    readTimer.start();
    parseBodFile(&parametersLocal,
		 &parametersWalkOnSpheres,
		 &parametersInteriorSampling,
		 &parametersResults,
		 &model);
    readTimer.stop();
  }
  
  Timer broadcastTimer;

  broadcastTimer.start();

  parametersWalkOnSpheres.mpiBroadcast(0);
  parametersInteriorSampling.mpiBroadcast(0);
  parametersVirial.mpiBroadcast(0);
  parametersResults.mpiBroadcast(0);
  parametersLocal.mpiBroadcast(0);
  model.mpiBroadcast(0);

  broadcastTimer.stop();

  if (parametersLocal.getXyzInputFileNameWasSet() !=
      parametersLocal.getMapInputFileNameWasSet()) {

    std::cerr << "Error: Both xyz and map file names must be given when using "
              << "trajectory mode" << std::endl;

    return 1;
  }

  if (parametersLocal.getXyzInputFileNameWasSet() &&
      !model.isEmpty()) {

    std::cerr << "Error: Cannot specify non-trajectory geometry when using "
	      << "trajectory mode" << std::endl;

    return 1;
  }

  if (false &&
      !(parametersInteriorSampling.getTotalNumSamplesWasSet() |
	parametersInteriorSampling.getMaxErrorVolumeWasSet())) {

    std::cerr << "Error: Must specify number of interior samples or volume "
	      << "error if computing form factors" << std::endl;

    return 1;
  }

  // Standard mode
  if (!parametersLocal.getXyzInputFileNameWasSet()) {
    if (parametersLocal.getPrintBenchmarks() && 
	parametersLocal.getMpiRank() == 0) {

      printRAM("loading input data",
	       &globalCsvItems);
    }

    perRunCsvItemsList.emplace_back();
    
    int runZenoStatus =
      runZeno(parametersLocal,
	      &parametersWalkOnSpheres,
	      &parametersInteriorSampling,
	      &parametersVirial,
	      &parametersResults,
	      readTimer.getTime(),
	      broadcastTimer.getTime(),
	      &model,
	      &perRunCsvItemsList.back());

    if (runZenoStatus != 0) {
      return runZenoStatus;
    }

    writeCsvOutputCols(globalCsvItems,
		       perRunCsvItemsList,
		       &csvOutputFile);
  }
  // Trajectory mode
  else {
    std::list<zeno::MixedModel<double> > snapshots;

    std::list<std::string> comments;
    
    if (parametersLocal.getMpiRank() == 0) {
      readTimer.start();
    
      std::unordered_map<std::string, double> atomIdToRadius;
    
      parseMapFile(parametersLocal,
		   &atomIdToRadius);
    
      parseXyzFile(parametersLocal,
		   atomIdToRadius,
		   &snapshots,
		   &comments);

      readTimer.stop();

      assert(snapshots.size() == comments.size());
    }

    broadcastTimer.start();

    broadcastSnapshots(&snapshots);

    // Only the rank 0 node prints the comments, but all nodes need to iterate
    // over the list, so don't broadcast the comments but resize the list
    comments.resize(snapshots.size());
    
    broadcastTimer.stop();

    if (parametersLocal.getPrintBenchmarks() && 
	parametersLocal.getMpiRank() == 0) {

      printRAM("loading input data",
	       &globalCsvItems);
    }

    auto commentsIterator = comments.begin();
    
    for (zeno::MixedModel<double> & snapshot: snapshots) {
      ParametersWalkOnSpheres
	snapshotParametersWalkOnSpheres(parametersWalkOnSpheres);
      ParametersInteriorSampling
	snapshotParametersInteriorSampling(parametersInteriorSampling);
      ParametersVirial
	snapshotParametersVirial(parametersVirial);
      // virial cannot be computed from snapshots
      snapshotParametersVirial.setSteps(0);
      ParametersResults
	snapshotParametersResults(parametersResults);

      perRunCsvItemsList.emplace_back();
      
      int runZenoStatus =
	runZeno(parametersLocal,
		&snapshotParametersWalkOnSpheres,
		&snapshotParametersInteriorSampling,
		&snapshotParametersVirial,
		&snapshotParametersResults,
		readTimer.getTime(),
		broadcastTimer.getTime(),
		&snapshot,
		&perRunCsvItemsList.back());

      if (runZenoStatus != 0) {
	return runZenoStatus;
      }

      printExactScalar("Comment", "comment", "",
		       *commentsIterator,
		       &perRunCsvItemsList.back());	

      ++ commentsIterator;
    }

    writeCsvOutputCols(globalCsvItems,
		       perRunCsvItemsList,
		       &csvOutputFile);

  }

  csvOutputFile.close();

#ifdef USE_MPI
  MPI_Finalize();
#endif
  
  if (mpiRank == 0) {
    printTime("End time: ");
  }
  
  return 0;
}

/// Parses the command line and stores the values of the parameters it
/// contains.  If a config file is specified on the command line, this will be
/// parsed as well.
///
/// Also computes default values for numThreads and seed if they were not set.
///
void 
parseCommandLine(int argc, char **argv,
		 ParametersWalkOnSpheres * parametersWalkOnSpheres,
		 ParametersInteriorSampling * parametersInteriorSampling,
		 ParametersVirial * parametersVirial,
		 ParametersResults * parametersResults,
		 ParametersLocal * parametersLocal) {
  
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
  
  parametersLocal->setInputFileName(args_info.input_file_arg);

  if (args_info.csv_output_file_given) {
    parametersLocal->setCsvOutputFileName(args_info.csv_output_file_arg);
  }
  
  parametersWalkOnSpheres->setMinTotalNumWalks(args_info.min_num_walks_arg);

  parametersInteriorSampling->setMinTotalNumSamples
    (args_info.min_num_interior_samples_arg);

  if (args_info.num_walks_given) {
    parametersWalkOnSpheres->setTotalNumWalks(args_info.num_walks_arg);
  }

  if (args_info.num_interior_samples_given) {
    parametersInteriorSampling->setTotalNumSamples
      (args_info.num_interior_samples_arg);
  }

  if (args_info.virial_steps_given) {
    parametersVirial->setSteps(args_info.virial_steps_arg);
  }

  if (args_info.max_rsd_capacitance_given) {
    parametersWalkOnSpheres->setMaxErrorCapacitance
      (args_info.max_rsd_capacitance_arg);
  }

  if (args_info.max_rsd_polarizability_given) {
    parametersWalkOnSpheres->setMaxErrorPolarizability
      (args_info.max_rsd_polarizability_arg);
  }

  if (args_info.max_rsd_volume_given) {
    parametersInteriorSampling->setMaxErrorVolume(args_info.max_rsd_volume_arg);
  }

  if (args_info.max_run_time_given) {
    parametersWalkOnSpheres->setMaxRunTime(args_info.max_run_time_arg);
    
    parametersInteriorSampling->setMaxRunTime(args_info.max_run_time_arg);
  }

  if (args_info.num_threads_given) {
    parametersWalkOnSpheres->setNumThreads(args_info.num_threads_arg);

    parametersInteriorSampling->setNumThreads(args_info.num_threads_arg);

    parametersVirial->setNumThreads(args_info.num_threads_arg);
  }
  else {
    int numThreads = std::thread::hardware_concurrency();
    
    if (numThreads == 0) {
      numThreads = 1;
    }
    
    parametersWalkOnSpheres->setNumThreads(numThreads);

    parametersInteriorSampling->setNumThreads(numThreads);

    parametersVirial->setNumThreads(numThreads);
  }

  if (args_info.seed_given) {
    parametersWalkOnSpheres->setSeed(args_info.seed_arg);

    parametersInteriorSampling->setSeed(args_info.seed_arg);

    parametersVirial->setSeed(args_info.seed_arg);
  }
  else {
    int seed{};
    
    int status = getRandomNumberFromOS(&seed);
      
    if (status != 0) {
      std::cerr << "Error randomly setting seed for random number generator"
		<< std::endl;
    }

    parametersWalkOnSpheres->setSeed(seed);

    parametersInteriorSampling->setSeed(seed);

    parametersVirial->setSeed(seed);
  }

  if (args_info.surface_points_file_given) {
    parametersLocal->setSurfacePointsFileName
      (args_info.surface_points_file_arg);

    parametersWalkOnSpheres->setSaveSurfacePoints(true);
  }

  if (args_info.interior_points_file_given) {
    parametersLocal->setInteriorPointsFileName
      (args_info.interior_points_file_arg);

    parametersInteriorSampling->setSaveInteriorPoints(true);
  }

  if (args_info.virial_coefficient_order_given) {
    parametersVirial->setOrder(args_info.virial_coefficient_order_arg);
  }

  parametersLocal->setPrintCounts(args_info.print_counts_given);
  parametersLocal->setPrintBenchmarks(args_info.print_benchmarks_given);

  free(parser_params);
}

/// Parses the bod file given as input and extracts sphere data, voxel data,
/// and parameters.
///
void
parseBodFile(ParametersLocal * parametersLocal,
	     ParametersWalkOnSpheres * parametersWalkOnSpheres,
	     ParametersInteriorSampling * parametersInteriorSampling,
	     ParametersResults * parametersResults,
	     MixedModel<double> * model) {

  std::string fileName = parametersLocal->getInputFileName();

  std::ifstream inputFile;

  inputFile.open(fileName, std::ifstream::in);

  if (!inputFile.is_open()) {
    std::cerr << "Error opening bod input file " << fileName << std::endl;
    exit(EXIT_FAILURE);
  }

  bod_parser::BodParser parser(parametersLocal,
			       inputFile,
			       parametersWalkOnSpheres,
			       parametersInteriorSampling,
			       parametersResults,
			       model);

  int parseResult = parser.parse();

  if (parseResult != 0) {
    std::cerr << "Error parsing bod input file " << fileName << std::endl;
    exit(EXIT_FAILURE);
  }
  
  inputFile.close();
}

/// Parses the map file given as input and extracts the atomIdToRadius map.
///
void
parseMapFile(ParametersLocal const & parametersLocal,
	     std::unordered_map<std::string, double> * atomIdToRadius) {
  
  std::string fileName = parametersLocal.getMapInputFileName();

  std::ifstream inputFile;

  inputFile.open(fileName, std::ifstream::in);

  if (!inputFile.is_open()) {
    std::cerr << "Error opening map input file " << fileName << std::endl;
    exit(EXIT_FAILURE);
  }

  map_parser::MapParser parser(inputFile,
			       atomIdToRadius);
  
  int parseResult = parser.parse();

  if (parseResult != 0) {
    std::cerr << "Error parsing map input file " << fileName << std::endl;
    exit(EXIT_FAILURE);
  }

  inputFile.close();
}

/// Parses the xyz file given as input and extracts the spheres representing
/// the atoms and any comment in each snapshot.
///
void
parseXyzFile(ParametersLocal const & parametersLocal,
	     std::unordered_map<std::string, double> const & atomIdToRadius,
	     std::list<zeno::MixedModel<double> > * snapshots,
	     std::list<std::string> * comments) {
  
  std::string fileName = parametersLocal.getXyzInputFileName();

  std::ifstream inputFile;

  inputFile.open(fileName, std::ifstream::in);

  if (!inputFile.is_open()) {
    std::cerr << "Error opening xyz input file " << fileName << std::endl;
    exit(EXIT_FAILURE);
  }


  xyz_parser::XyzParser parser(atomIdToRadius,
			       inputFile,
			       snapshots,
			       comments);

  int parseResult = parser.parse();

  if (parseResult != 0) {
    std::cerr << "Error parsing xyz input file " << fileName << std::endl;
    exit(EXIT_FAILURE);
  }
  
  inputFile.close();
}

/// If MPI rank is 0, broadcasts the list of snapshots to all other MPI nodes.
/// If MPI rank is not 0, fills the list of snapshots with the snapshots
/// received from the rank 0 node.
///
void
broadcastSnapshots(std::list<zeno::MixedModel<double> > * snapshots) {
#ifdef USE_MPI
  unsigned long long int numSnapshots = snapshots->size();

  MPI_Bcast(&numSnapshots, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);

  snapshots->resize(numSnapshots);

  for (zeno::MixedModel<double> & snapshot: *snapshots) {
    snapshot.mpiBroadcast(0);
  }
#endif
}

/// Runs the Zeno computations on the given model and outputs the results.
///
/// Returns 0 on success, and 1 otherwise.
///
int
runZeno(ParametersLocal const & parametersLocal,
	ParametersWalkOnSpheres * parametersWalkOnSpheres,
	ParametersInteriorSampling * parametersInteriorSampling,
	ParametersVirial * parametersVirial,
	ParametersResults * parametersResults,
	double readTime,
	double broadcastTime,
	MixedModel<double> * model,
	CsvItems * csvItems) {
  
  Zeno zeno(model);

  if (parametersLocal.getPrintBenchmarks() && 
      parametersLocal.getMpiRank() == 0) {

    printRAM("building spatial data structure",
	     csvItems);
  }

  Zeno::Status doWalkOnSpheresStatus =
    zeno.doWalkOnSpheres(parametersWalkOnSpheres,
			 parametersResults);

  if (doWalkOnSpheresStatus != Zeno::Status::Success) {
    if (doWalkOnSpheresStatus == Zeno::Status::EmptyModel) {
      std::cerr << "Error: no geometry loaded" << std::endl;
    }
    
    return 1;
  }

  if (parametersLocal.getPrintBenchmarks() && 
      parametersLocal.getMpiRank() == 0) {

    printRAM("walk on spheres",
	     csvItems);
  }

  Zeno::Status doInteriorSamplingStatus =
    zeno.doInteriorSampling(parametersInteriorSampling,
			    parametersResults);

  if (doInteriorSamplingStatus != Zeno::Status::Success) {
    if (doInteriorSamplingStatus == Zeno::Status::EmptyModel) {
      std::cerr << "Error: no geometry loaded" << std::endl;
    }
    
    return 1;
  }

  if (parametersLocal.getPrintBenchmarks() && 
      parametersLocal.getMpiRank() == 0) {

    printRAM("interior samples",
	     csvItems);
  }
  
  Zeno::Status doVirialSamplingStatus =
    zeno.doVirialSampling(parametersVirial,
			  parametersResults);

  if (doVirialSamplingStatus != Zeno::Status::Success) {
    if (doVirialSamplingStatus == Zeno::Status::EmptyModel) {
      std::cerr << "Error: no geometry loaded" << std::endl;
    }
    
    return 1;
  }

  if (parametersLocal.getPrintBenchmarks() && 
      parametersLocal.getMpiRank() == 0) {

    printRAM("interior samples",
	     csvItems);
  }
  
  if (parametersWalkOnSpheres->getMaxRunTimeWasSet() &&
      (zeno.getTotalTime() > parametersWalkOnSpheres->getMaxRunTime())) {

    std::cerr << std::endl
              << "*** Warning ***" << std::endl
	      << "Max run time was reached.  Not all requested computations "
	      << "may have been performed." << std::endl;
  }

  Results results;

  zeno.getResults(parametersResults,
		  &results);
  
  printOutput(results,
	      *parametersWalkOnSpheres,
	      *parametersInteriorSampling,
	      *parametersVirial,
	      *parametersResults,
	      parametersLocal,
	      zeno.getInitializeTime(),
	      readTime,
	      broadcastTime,
	      zeno.getPreprocessTime(),
	      zeno.getWalkOnSpheresTime(),
	      zeno.getWalkOnSpheresReductionTime(),
	      zeno.getInteriorSamplingTime(),
	      zeno.getInteriorSamplingReductionTime(),
	      zeno.getVirialTime(),
	      zeno.getVirialReductionTime(),
	      zeno.getTotalTime(),
	      csvItems);
  
  savePointFiles(&zeno,
		 parametersLocal);
  
  return 0;
}

/// Prints parameters, results, and (optionally) detailed running time
/// benchmarks on MPI process 0.
///
void
printOutput(Results const & results,
	    ParametersWalkOnSpheres const & parametersWalkOnSpheres,
	    ParametersInteriorSampling const & parametersInteriorSampling,
	    ParametersVirial const & parametersVirial,
	    ParametersResults const & parametersResults,
	    ParametersLocal const & parametersLocal,
	    double initializeTime,
	    double readTime,
	    double broadcastTime,
	    double preprocessTime,
	    double walkTime,
	    double reduceTime,
	    double sampleTime,
	    double volumeReduceTime,
	    double virialTime,
	    double virialReduceTime,
	    double totalZenoTime,
	    CsvItems * csvItems) {
  
  if (parametersLocal.getMpiRank() == 0) {    
    std::cout << std::endl
	      << "Parameters" << std::endl
	      << "----------" << std::endl
	      << std::endl;

    printParameters(parametersWalkOnSpheres,
		    parametersInteriorSampling,
		    parametersVirial,
		    parametersResults,
		    parametersLocal,
		    csvItems);

    std::cout << std::endl
	      << "Results" << std::endl
	      << "-------" << std::endl
	      << std::endl;

    printResults(results,
		 csvItems);

    if (parametersLocal.getPrintCounts()) {
      std::cout << "Counts:" << std::endl
		<< std::endl;

      if (results.resultsZenoCompiled) {
	printScalar(results.t,
		    csvItems);

	printVector3(results.u,
		     csvItems);

	printMatrix3x3(results.v,
		       csvItems);

	printMatrix3x3(results.w,
		       csvItems);
      }

      if (results.resultsInteriorCompiled) {
	printScalar(results.numInteriorHits,
		    csvItems);
      }
    } 

    if (parametersLocal.getPrintBenchmarks()) {
      printExactScalar("Initialize     ",
		       "initialize_time",
		       "s",
		       initializeTime,	      
		       csvItems);

      printExactScalar("Read           ",
		       "read_time",
		       "s",
		       readTime,
		       csvItems);

      printExactScalar("Broadcast      ",
		       "broadcast_time",
		       "s",
		       broadcastTime,
		       csvItems);
      
      printExactScalar("Preprocess     ",
		       "preprocess_time",
		       "s",
		       preprocessTime,
		       csvItems);

      printExactScalar("Exterior Walk  ",
		       "exterior_walk_time",
		       "s",
		       walkTime,
		       csvItems);

      printExactScalar("Exterior Reduce",
		       "exterior_reduce_time",
		       "s",
		       reduceTime,
		       csvItems);

      printExactScalar("Volume Sample  ",
		       "volume_sample_time",
		       "s",
		       sampleTime,
		       csvItems);

      printExactScalar("Volume Reduce  ",
		       "volume_reduce_time",
		       "s",
		       volumeReduceTime,
		       csvItems);

      printExactScalar("Virial Sample  ",
		       "virial_sample_time",
		       "s",
		       virialTime,
		       csvItems);

      printExactScalar("Virial Reduce  ",
		       "virial_reduce_time",
		       "s",
		       virialReduceTime,
		       csvItems);

      printExactScalar("Total Time     ",
		       "total_time",
		       "s",
		       totalZenoTime,
		       csvItems);
      
      std::cout << std::endl;
    }
  }
}

/// Prints the parameters.  Most parameters are not printed if they have not
/// been set.
/// 
void
printParameters(ParametersWalkOnSpheres const & parametersWalkOnSpheres,
	        ParametersInteriorSampling const & parametersInteriorSampling,
	        ParametersVirial const & parametersVirial,
	        ParametersResults const & parametersResults,
		ParametersLocal const & parametersLocal,
		CsvItems * csvItems) {
  
  printExactScalar("Input file", "input_file", "",
		   parametersLocal.getInputFileName(),
		   csvItems);

  printExactScalar("Number of nodes", "num_nodes", "",
		   parametersLocal.getMpiSize(),
		   csvItems);

  printExactScalar("Number of threads", "num_threads", "",
		   parametersWalkOnSpheres.getNumThreads(),
		   csvItems);

  printExactScalar("Random number seed", "random_number_seed", "",
		   parametersWalkOnSpheres.getSeed(),
		   csvItems);

  if (parametersWalkOnSpheres.getTotalNumWalksWasSet()) {
    printExactScalar("Number of walks requested", "num_walks_requested", "",
		     parametersWalkOnSpheres.getTotalNumWalks(),
		     csvItems);
  }

  if (parametersInteriorSampling.getTotalNumSamplesWasSet()) {
    printExactScalar("Number of interior samples requested",
		     "num_interior_samples_requested", "",
		     parametersInteriorSampling.getTotalNumSamples(),
		     csvItems);
  }

  if (parametersWalkOnSpheres.getMaxErrorCapacitanceWasSet()) {
    printExactScalar("Max error in capacitance requested",
		     "max_error_capacitance_requested", "%",
	             parametersWalkOnSpheres.getMaxErrorCapacitance(),
		     csvItems);
  }

  if (parametersWalkOnSpheres.getMaxErrorPolarizabilityWasSet()) {
    printExactScalar("Max error in mean polarizability requested",
		     "max_error_mean_polarizability_requested", "%",
	             parametersWalkOnSpheres.getMaxErrorPolarizability(),
		     csvItems);
  }

  if (parametersInteriorSampling.getMaxErrorVolumeWasSet()) {
    printExactScalar("Max error in volume requested",
		     "max_error_volume_requested", "%",
	             parametersInteriorSampling.getMaxErrorVolume(),
		     csvItems);
  }

  if (parametersWalkOnSpheres.getMaxRunTimeWasSet()) {
    printExactScalar("Max run time", "max_run_time", "s",
		     parametersWalkOnSpheres.getMaxRunTime(),
		     csvItems);
  }

  if (parametersWalkOnSpheres.getSkinThicknessWasSet()) {
    printExactScalar("Skin thickness", "skin_thickness", "",
		     parametersWalkOnSpheres.getSkinThickness(),
		     csvItems);
  }

  if (parametersWalkOnSpheres.getLaunchCenterWasSet()) {
    printExactVector3("Launch center", "launch_center", "",
		      parametersWalkOnSpheres.getLaunchCenter(),
		      csvItems);
  }

  if (parametersWalkOnSpheres.getLaunchRadiusWasSet()) {
    printExactScalar("Launch radius", "launch_radius", "",
		     parametersWalkOnSpheres.getLaunchRadius(),
		     csvItems);
  }

  if (parametersVirial.getOrderWasSet()) {
    printExactScalar("Virial coefficient order", "virial_order", "",
		     parametersVirial.getOrder(),
		     csvItems);
  }

  if (parametersVirial.getStepsWasSet()) {
    printExactScalar("Virial steps", "virial_steps", "",
		     parametersVirial.getSteps(),
		     csvItems);
  }



  if (parametersResults.getLengthScaleWasSet()) {
    printExactScalar("Length scale", "length_scale",
		     Units::getName(parametersResults.getLengthScaleUnit()),
		     parametersResults.getLengthScaleNumber(),
		     csvItems);
  }

  if (parametersResults.getTemperatureWasSet()) {
    printExactScalar("Temperature", "temperature",
		     Units::getName(parametersResults.getTemperatureUnit()),
		     parametersResults.getTemperatureNumber(),
		     csvItems);
  }

  if (parametersResults.getMassWasSet()) {
    printExactScalar("Mass", "mass",
		     Units::getName(parametersResults.getMassUnit()),
		     parametersResults.getMassNumber(),
		     csvItems);
  }

  if (parametersResults.getSolventViscosityWasSet()) {
    printExactScalar("Solvent viscosity", "solvent_viscosity",
		     Units::getName
		     (parametersResults.getSolventViscosityUnit()),
		     parametersResults.getSolventViscosityNumber(),
		     csvItems);
  }

  if (parametersResults.getBuoyancyFactorWasSet()) {
    printExactScalar("Buoyancy factor", "buoyancy_factor", "",
		     parametersResults.getBuoyancyFactor(),
		     csvItems);
  }
}

/// Prints the results.  Results that were not computed are not printed.
///
void
printResults(Results const & results,
	     CsvItems * csvItems) {

  if (results.resultsZenoCompiled) {
    printExactScalar(results.numWalks.prettyName,
		     results.numWalks.csvName,
		     results.numWalks.unit,
		     results.numWalks.value,
		     csvItems);

    std::cout << std::endl;
  }

  if (results.resultsInteriorCompiled) {
    printExactScalar(results.numInteriorSamples.prettyName,
		     results.numInteriorSamples.csvName,
		     results.numInteriorSamples.unit,
		     results.numInteriorSamples.value,
		     csvItems);

    std::cout << std::endl;
  }
  
  if (results.resultsVirialCompiled) {
    printExactScalar(results.refFrac.prettyName,
		     results.refFrac.csvName,
		     results.refFrac.unit,
		     results.refFrac.value,
		     csvItems);

    printScalar(results.virialCoefficient,
		csvItems);

    std::cout << std::endl;
  }

  if (results.resultsZenoCompiled) {

    printScalar(results.capacitance,
		csvItems);

    printMatrix3x3(results.polarizabilityTensor,
		   csvItems);

    printVector3(results.polarizabilityEigenvalues,
		 csvItems);

    printScalar(results.meanPolarizability,
		csvItems);

    printScalar(results.hydrodynamicRadius,
		csvItems);

    printScalar(results.q_eta,
		csvItems);

    printScalar(results.viscometricRadius,
		csvItems);

    if (results.intrinsicViscosityConventionalComputed) {
      
      printScalar(results.intrinsicViscosityConventional,
		  csvItems);
    }

    if (results.frictionCoefficientComputed) {

      printScalar(results.frictionCoefficient,
		  csvItems);
    }

    if (results.diffusionCoefficientComputed) {

      printScalar(results.diffusionCoefficient,
		  csvItems);
    }

    if (results.sedimentationCoefficientComputed) {

      printScalar(results.sedimentationCoefficient,
		  csvItems);
    }
  }

  if (results.resultsInteriorCompiled) {

    printScalar(results.volume,
		csvItems);

    printScalar(results.capacitanceOfASphere,
		csvItems);

    printMatrix3x3(results.gyrationTensor,
		   csvItems);

    printVector3(results.gyrationEigenvalues,
		 csvItems);
  }

  if (results.resultsZenoCompiled &&
      results.resultsInteriorCompiled) {

    printScalar(results.intrinsicConductivity,
		csvItems);

    printScalar(results.intrinsicViscosity,
		csvItems);
  }

  // if (results.formResultsCompiled) {
  //   if (!csvFormat) {
  //     *out << "Form factor: " << std::endl;

  //     for (unsigned int factorNum = 0; 
  // 	   factorNum < numFormFactors;
  // 	   ++factorNum) {

  // 	*out << std::scientific
  // 	     << formFactorQs.at(factorNum) << " ("
  // 	     << Units::getName(parameters->getLengthScaleUnit())
  // 	     << "^-1): "
  // 	     << formFactors.at(factorNum) << std::endl;
  //     }

  //     *out << std::endl;
  //   }
  // }
}

/// Prints a scalar that does not have uncertainty
///
template <typename T>
void
printExactScalar(std::string const & prettyName,
		 std::string const & csvName,
		 std::string const & units,
		 T property,
		 CsvItems * csvItems) {

  std::cout << prettyName;

  if (!units.empty()) {
    std::cout << " (" << units << ")";
  }

  std::cout << std::fixed << ": " << property << std::endl;

  if (!units.empty()) {
    csvItems->at(0).push_back(csvName);
    csvItems->at(1).push_back("units");
    csvItems->at(2).push_back(units);
  }

  csvItems->at(0).push_back(csvName);
  csvItems->at(1).push_back("value");
  csvItems->at(2).push_back(to_string_scientific(property));
}

/// Prints a Vector3 that does not have uncertainty
///
template <typename T>
void
printExactVector3(std::string const & prettyName,
		  std::string const & csvName,
		  std::string const & units,
		  Vector3<T> const & property,
		  CsvItems * csvItems) {

  std::cout << prettyName;

  if (!units.empty()) {
    std::cout << " (" << units << ")";
  }

  std::cout << std::fixed << ": " << property << std::endl;
  
  for (int i = 0; i < 3; ++i) {
    if (!units.empty()) {
      csvItems->at(0).push_back(csvName + "[" + std::to_string(i) + "]");
      csvItems->at(1).push_back("units");
      csvItems->at(2).push_back(units);
    }

    csvItems->at(0).push_back(csvName + "[" + std::to_string(i) + "]");
    csvItems->at(1).push_back("value");
    csvItems->at(2).push_back(to_string_scientific(property.get(i)));
  }
}

/// Prints a scalar with uncertainty
///
void
printScalar(Result<Uncertain<double> > const & result,
	    CsvItems * csvItems) {

  std::cout << std::scientific
	    << result.prettyName << " (" << result.unit << "): "
	    << result.value << std::endl
	    << std::endl;

  csvItems->at(0).push_back(result.csvName);
  csvItems->at(1).push_back("units");
  csvItems->at(2).push_back(result.unit);

  csvItems->at(0).push_back(result.csvName);
  csvItems->at(1).push_back("value");
  csvItems->at(2).push_back(to_string_scientific(result.value.getMean()));

  csvItems->at(0).push_back(result.csvName);
  csvItems->at(1).push_back("std_dev");
  csvItems->at(2).push_back(to_string_scientific(result.value.getStdDev()));
}

/// Prints a Vector3 with uncertainty
///
void
printVector3(Result<Vector3<Uncertain<double> > > const & result,
	     CsvItems * csvItems) {

  std::cout << std::scientific
	    << result.prettyName << " (" << result.unit << "): "
	    << std::endl
	    << result.value << std::endl
	    << std::endl;

  for (int i = 0; i < 3; ++i) {
    csvItems->at(0).push_back(result.csvName + "[" + std::to_string(i) + "]");
    csvItems->at(1).push_back("units");
    csvItems->at(2).push_back(result.unit);

    csvItems->at(0).push_back(result.csvName + "[" + std::to_string(i) + "]");
    csvItems->at(1).push_back("value");
    csvItems->at(2).push_back(to_string_scientific(result.value.get(i).getMean()));

    csvItems->at(0).push_back(result.csvName + "[" + std::to_string(i) + "]");
    csvItems->at(1).push_back("std_dev");
    csvItems->at(2).push_back(to_string_scientific(result.value.get(i).getStdDev()));
  }
}

/// Prints a Matrix3x3 with uncertainty
///
void
printMatrix3x3(Result<Matrix3x3<Uncertain<double> > > const & result,
	       CsvItems * csvItems) {

  std::cout << std::scientific
	    << result.prettyName << " (" << result.unit << "): "
	    << std::endl
	    << result.value << std::endl
	    << std::endl;

  for (int row = 0; row < 3; ++row) {
    for (int col = 0; col < 3; ++col) {
      csvItems->at(0).push_back(result.csvName + "[" + std::to_string(row) + "][" + std::to_string(col) + "]");
      csvItems->at(1).push_back("units");
      csvItems->at(2).push_back(result.unit);

      csvItems->at(0).push_back(result.csvName + "[" + std::to_string(row) + "][" + std::to_string(col) + "]");
      csvItems->at(1).push_back("value");
      csvItems->at(2).push_back(to_string_scientific(result.value.get(row, col).getMean()));

      csvItems->at(0).push_back(result.csvName + "[" + std::to_string(row) + "][" + std::to_string(col) + "]");
      csvItems->at(1).push_back("std_dev");
      csvItems->at(2).push_back(to_string_scientific(result.value.get(row, col).getStdDev()));
    }
  }
}

/// Writes Walk-on-Spheres and Interior Sampling hit points to disk.
/// 
void
savePointFiles(Zeno * zeno,
	       ParametersLocal const & parametersLocal) {

  if (!parametersLocal.getSurfacePointsFileName().empty()) {
    std::vector<Vector3<double> > const * points = nullptr;
    std::vector<Vector3<char> > const * charges = nullptr;
    
    zeno->getWalkOnSpheresHitPoints(&points, &charges);
    
    if (points == nullptr ||
	charges == nullptr) {
      
      std::cerr << "*** Warning ***" << std::endl
		<< "A surface points file was requested but walks were not "
		<< "performed.  Surface points file will not be written."
		<< std::endl
		<< std::endl;
    }
    else {
      writePoints(parametersLocal.getSurfacePointsFileName(), 
		  points, 
		  charges);
    }
  }

  if (!parametersLocal.getInteriorPointsFileName().empty()) {
    std::vector<Vector3<double> > const * points = nullptr;

    zeno->getInteriorSamplingHitPoints(&points);
    
    if (points == nullptr) {
      
      std::cerr << "*** Warning ***" << std::endl
		<< "An interior points file was requested but interior "
		<< "samples were not performed.  Interior points file will "
		<< "not be written."
		<< std::endl
		<< std::endl;
    }
    else {
      writePoints(parametersLocal.getInteriorPointsFileName(), 
		  points, 
		  nullptr);
    }
  }
}

/// Prints the current date and time prefixed by the given label.
///
void
printTime(std::string const & label) {
  const int bufferSize = 256;

  time_t rawtime;
  struct tm * timeinfo;
  char buffer[bufferSize];
    
  time(&rawtime);
  timeinfo = localtime(&rawtime);
    
  strftime(buffer, bufferSize, "%F %T", timeinfo);
  
  std::cout << label << buffer << std::endl;
}

/// Prints the current RAM used by the process prefixed by the given label.
///
void
printRAM(std::string const & label,
	 CsvItems * csvItems) {
  
  std::ifstream statusFile("/proc/self/status");

  std::string line;

  size_t pos = std::string::npos;

  while (statusFile.good() &&
	 pos == std::string::npos) {

    std::getline(statusFile, line);

    pos = line.find("VmRSS:");
  }

  statusFile.close();

  std::string ram;
  std::string units;
  
  if (pos != std::string::npos) {
    size_t unitsPos = line.find_last_of(" ") + 1;

    size_t ramPos = line.find_last_of("\t :", unitsPos - 2) + 1;
      
    ram = line.substr(ramPos, unitsPos - ramPos);
    
    units = line.substr(unitsPos);
  }
  else {
    ram = "unknown";
    units = "";
  }

  std::string prettyName = "RAM after " + label;
  
  std::string csvName = label + "_ram";
  std::replace(csvName.begin(), csvName.end(), ' ', '_');
  
  printExactScalar(prettyName,
		   csvName,
		   units,
		   ram,	       
		   csvItems);
}

/// Writes the given CSV items into the given CSV output file.
/// The header/type/data tuples are written as rows down a column.
///
void
writeCsvOutputCols(CsvItems const & globalCsvItems,
		   std::list<CsvItems> const & perRunCsvItemsList,
		   std::ofstream * csvOutputFile) {

  assert(globalCsvItems.at(0).size() == globalCsvItems.at(1).size());
  assert(globalCsvItems.at(0).size() == globalCsvItems.at(2).size());

  // write global CSV items
  for (size_t row = 0;
       row < globalCsvItems.at(0).size();
       ++row) {

    // write header
    *csvOutputFile << globalCsvItems.at(0).at(row) << ","
		   << globalCsvItems.at(1).at(row);

    // write data
    for (size_t i = 0;
	 i < perRunCsvItemsList.size();
	 ++i) {
      
      *csvOutputFile << "," << globalCsvItems.at(2).at(row);
    }

    *csvOutputFile << std::endl;
  }
  
  // write per-run CSV items
  for (size_t row = 0;
       row < perRunCsvItemsList.back().at(0).size();
       ++row) {

    // write header
    *csvOutputFile << perRunCsvItemsList.back().at(0).at(row) << ","
		   << perRunCsvItemsList.back().at(1).at(row);

    // write data
    for (CsvItems const & perRunCsvItems: perRunCsvItemsList) {
      assert(perRunCsvItems.at(0).size() == perRunCsvItems.at(1).size());
      assert(perRunCsvItems.at(0).size() == perRunCsvItems.at(2).size());

      assert(perRunCsvItems.at(0).size() == perRunCsvItemsList.back().at(0).size());
    
      *csvOutputFile << "," << perRunCsvItems.at(2).at(row);
    }

    *csvOutputFile << std::endl;
  }
}

/// Writes the given CSV items into the given CSV output file.
/// The header/type/data tuples are written as columns across a row.
///
void
writeCsvOutputRows(CsvItems const & globalCsvItems,
		   std::list<CsvItems> const & perRunCsvItemsList,
		   std::ofstream * csvOutputFile) {

  bool writeHeader = true;
  
  for (CsvItems const & perRunCsvItems: perRunCsvItemsList) {
    
    size_t row = 0;

    if (!writeHeader) {
      row = 2;
    }

    for (; row < 3; ++row) {
      bool writeComma = false;
    
      assert(globalCsvItems.at(0).size() == globalCsvItems.at(1).size());
      assert(globalCsvItems.at(0).size() == globalCsvItems.at(2).size());
    
      for (size_t i = 0; i < globalCsvItems.at(row).size(); ++i) {
	if (writeComma) {
	  *csvOutputFile << ",";
	}
    
	*csvOutputFile << globalCsvItems.at(row).at(i);

	writeComma = true;
      }

      assert(perRunCsvItems.at(0).size() == perRunCsvItems.at(1).size());
      assert(perRunCsvItems.at(0).size() == perRunCsvItems.at(2).size());

      assert(perRunCsvItems.at(0).size() == perRunCsvItemsList.back().at(0).size());
      
      for (size_t i = 0; i < perRunCsvItems.at(row).size(); ++i) {
	if (writeComma) {
	  *csvOutputFile << ",";
	}
    
	*csvOutputFile << perRunCsvItems.at(row).at(i);

	writeComma = true;
      }

      *csvOutputFile << std::endl;
    }

    writeHeader = false;
  }
}

/// Writes the given set of points and (if not NULL) charges to disk.
///
void
writePoints(std::string const & fileName, 
	    std::vector<Vector3<double> > const * points,
	    std::vector<Vector3<char> > const * charges) {

  std::ofstream outputFile;

  outputFile.open(fileName, std::ofstream::out);

  if (!outputFile.is_open()) {

    std::cerr << "Error opening output file " << fileName << std::endl;
    exit(1);
  }

  for (unsigned int i = 0; i < points->size(); i++) {

    if (charges != nullptr) {
      outputFile << charges->at(i).get(0)
		 << charges->at(i).get(1)
		 << charges->at(i).get(2);
    }

    outputFile << std::setw(16) << points->at(i).get(0)
	       << std::setw(16) << points->at(i).get(1)
	       << std::setw(16) << points->at(i).get(2) << std::endl;
  }

  outputFile.close();
}

/// Sets the variable pointed to by "randomNumber" to a random number obtained
/// from the OS.
///
/// Returns 0 on success, and 1 otherwise.
///
template <class T>
int
getRandomNumberFromOS(T * randomNumber) {
  int status = 0;
    
  std::ifstream urandom("/dev/urandom");

  urandom.read((char *)randomNumber, sizeof(*randomNumber));

  if (!urandom.good()) {
    status = 1;
  }

  urandom.close();

  return status;
}

/// Converts a value to a string in the "scientific" format
///
template <typename T>
std::string
to_string_scientific(T const & value)
{
  std::ostringstream out;
  out << std::scientific << value;
  return out.str();
}
