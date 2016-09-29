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
// Date:    Mon Feb 24 11:25:18 2014 EDT
//
// Time-stamp: <2016-09-28 12:08:25 dcj>
//
// ================================================================

#ifdef USE_MPI
#include <mpi.h>
#endif

#include <thread>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <ctime>

#include "Parameters.h"
#include "Units.h"
#include "Uncertain.h"

#include "Parser.h"

#include "ResultsZeno.h"
#include "ResultsInterior.h"
#include "ResultsCompiler.h"

#include "Geometry/Sphere.h"
#include "Geometry/Spheres.h"
#include "Geometry/Vector3.h"

#include "NearestSurfacePoint/PointFromSphereCenters.h"

#include "InsideOutside/InOutSphereCenters.h"

#include "BoundingSphere/BoundingSphereAABB.h"

#include "SpherePoint/RandomSpherePointMarsaglia.h"
#include "SpherePoint/BiasedSpherePointDirect.h"
#include "SpherePoint/RandomBallPointRejection.h"

#include "Walker/WalkerExterior.h"
#include "Walker/SamplerInterior.h"

#include "Timer.h"

#include "SphereCenterModel/NanoFLANNSort.h"

#include "RandomNumber/SPRNG.h"

using SpheresModel = NanoFLANNSort;

using SpheresNearestSurfacePointFinder = PointFromSphereCenters<SpheresModel>;
using SpheresInsideOutsideTester       = InOutSphereCenters<SpheresModel>;

using RandomNumberGenerator = SPRNG;

using BoundingSphereGenerator = 
  BoundingSphereAABB<double>;

using RandomSpherePointGenerator = 
  RandomSpherePointMarsaglia<double, 
			     RandomNumberGenerator>;

using BiasedSpherePointGenerator = 
  BiasedSpherePointDirect<double, 
			  RandomNumberGenerator, 
			  RandomSpherePointGenerator>;

using RandomBallPointGenerator =
  RandomBallPointRejection<double,
			   RandomNumberGenerator>;

// ================================================================

int
getInput(int argc, char **argv,
	 Parameters * parameters,
	 long long * numWalksInProcess,
	 long long * numSamplesInProcess,
	 Spheres<double> * spheres,
	 bool * spheresLoaded,
	 double * initializeTime,
	 double * readTime,
	 double * broadcastTime);

int
preprocessWalkOnSpheres(bool spheresLoaded,
			Parameters const & parameters,
			Spheres<double> const & spheres,
			Sphere<double> * boundingSphere,
			SpheresModel * spheresModel,
			SpheresNearestSurfacePointFinder * *
			spheresNearestSurfacePointFinder,
			double * preprocessTime);

void setupRNGs(Parameters const & parameters,
	       std::vector<RandomNumberGenerator> * threadRNGs);

int
getWalkOnSpheresResults(long long numWalksInProcess,			
			Parameters const & parameters,
			Sphere<double> const & boundingSphere,
			SpheresNearestSurfacePointFinder const * 
			spheresNearestSurfacePointFinder,			
			std::vector<RandomNumberGenerator> * threadRNGs,
			ResultsZeno * * resultsZeno,
			double * walkTime,
			double * reduceTime);

int
preprocessInterior(bool spheresLoaded,
		   Parameters const & parameters,
		   SpheresModel const & spheresModel,
		   SpheresInsideOutsideTester * *
		   spheresInsideOutsideTester,
		   double * surfacePreprocessTime);

int
getInteriorResults(long long numSamplesInProcess,			
		   Parameters const & parameters,
		   Sphere<double> const & boundingSphere,
		   SpheresInsideOutsideTester const * 
		   spheresInsideOutsideTester,			
		   std::vector<RandomNumberGenerator> * threadRNGs,
		   ResultsInterior * * resultsInterior,
		   double * sampleTime,
		   double * volumeReduceTime);

long long 
computeNumInProcess(int mpiSize, int mpiRank,
		    long long totalNumSamples);

void
parseBodFile(Parameters * parameters, 
	     Spheres<double> * spheres);

void
getBodData(Parameters * parameters,
	   Spheres<double> * spheres,
	   double * readTime,
	   double * broadcastTime);

void
computeDefaultParameters(Parameters * parameters,
			 Sphere<double> * boundingSphere);

long long 
estimateTotalNum(double requestedError,
		 long long numSoFar,
		 Uncertain<double> const & currentValue);

void 
doWalkOnSpheresSelector(Parameters const & parameters,
			long long numWalksInProcess,
			Sphere<double> const & boundingSphere, 
			SpheresNearestSurfacePointFinder const * 
			spheresNearestSurfacePointFinder,
			std::vector<RandomNumberGenerator> * threadRNGs,
			ResultsZeno * resultsZeno,
			double * walkTime);

template <class NearestSurfacePointFinder>
void
doWalkOnSpheres(int numThreads,
		long long numWalksInProcess,
		Sphere<double> const & boundingSphere, 
		NearestSurfacePointFinder const & 
		nearestSurfacePointFinder,
		double fracErrorBound,
		double shellThickness,
		std::vector<RandomNumberGenerator> * threadRNGs,
		ResultsZeno * resultsZeno,
		double * walkTime);

template <class NearestSurfacePointFinder>
void
doWalkOnSpheresThread(Sphere<double> const & boundingSphere, 
		      NearestSurfacePointFinder const & 
		      nearestSurfacePointFinder,
		      int threadNum,
		      double fracErrorBound,
		      double shellThickness,
		      long long numWalks,
		      RandomNumberGenerator * randomNumberGenerator,
		      ResultsZeno * resultsZeno);

void
doInteriorSamplingSelector(Parameters const & parameters,
			   long long numSamplesInProcess,
			   Sphere<double> const & boundingSphere,
			   SpheresInsideOutsideTester const *  
			   spheresInsideOutsideTester,
			   std::vector<RandomNumberGenerator> * threadRNGs,
			   ResultsInterior * resultsInterior,
			   double * sampleTime);

template <class InsideOutsideTester>
void
doInteriorSampling(int numThreads,
		   long long numSamplesInProcess,
		   Sphere<double> const & boundingSphere, 
		   InsideOutsideTester const & insideOutsideTester,
		   double fracErrorBound,
		   std::vector<RandomNumberGenerator> * threadRNGs,
		   ResultsInterior * resultsInterior,
		   double * sampleTime);

template <class InsideOutsideTester>
void
doInteriorSamplingThread(Sphere<double> const & boundingSphere, 
			 InsideOutsideTester const & insideOutsideTester,
			 int threadNum,
			 double fracErrorBound,
			 long long numSamples,
			 RandomNumberGenerator * randomNumberGenerators,
			 ResultsInterior * resultsInterior);

void
printOutput(Sphere<double> const & boundingSphere,
	    ResultsInterior const * resultsInterior,
	    ResultsZeno const * resultsZeno,
	    Parameters const & parameters, 
	    double initializeTime,
	    double readTime,
	    double broadcastTime,
	    double preprocessTime,
	    double walkTime,
	    double reduceTime,
	    double surfacePreprocessTime,
	    double sampleTime,
	    double volumeReduceTime);

void
savePointFiles(ResultsInterior & resultsInterior,
	       ResultsZeno & resultsZeno,
	       Parameters const & parameters);

void
printTime(std::string const & label);

void
printRAM(std::string const & label);

void
writePoints(std::string const & fileName, 
	    std::vector<Vector3<double> > const * points,
	    std::vector<Vector3<char> > const * charges);

// ================================================================

int main(int argc, char **argv) {

  Timer totalTimer;

  totalTimer.start();

  Timer initializeTimer;

  initializeTimer.start();

  int mpiSize = 1, mpiRank = 0;

#ifdef USE_MPI
  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
#endif

  if (mpiRank == 0) {
    printTime("Start time: ");
  }

  Parameters parameters;

  parameters.setMpiSize(mpiSize);
  parameters.setMpiRank(mpiRank);

  long long numWalksInProcess   = 0;
  long long numSamplesInProcess = 0;

  double initializeTime = 0;
  double readTime       = 0;
  double broadcastTime  = 0;

  Spheres<double> spheres;

  bool spheresLoaded = false;

  initializeTimer.stop();

  int getInputSuccess = 
    getInput(argc, argv,
	     &parameters,
	     &numWalksInProcess,
	     &numSamplesInProcess,
	     &spheres,
	     &spheresLoaded,
	     &initializeTime,
	     &readTime,
	     &broadcastTime);

  if (getInputSuccess != 0) {
    return getInputSuccess;
  }

  initializeTimer.start();

  SpheresModel spheresModel;

  SpheresNearestSurfacePointFinder * spheresNearestSurfacePointFinder = NULL;

  Sphere<double> boundingSphere;

  double preprocessTime = 0;

  initializeTimer.stop();

  int preprocessWalkOnSpheresSuccess = 
    preprocessWalkOnSpheres(spheresLoaded,
			    parameters,
			    spheres,
			    &boundingSphere,
			    &spheresModel,
			    &spheresNearestSurfacePointFinder,
			    &preprocessTime);

  if (preprocessWalkOnSpheresSuccess != 0) {
    return preprocessWalkOnSpheresSuccess;
  }

  initializeTimer.start();

  computeDefaultParameters(&parameters, &boundingSphere);

  std::vector<RandomNumberGenerator> threadRNGs;

  setupRNGs(parameters,
	    &threadRNGs);

  ResultsZeno * resultsZeno = NULL;

  double walkTime   = 0;
  double reduceTime = 0;

  initializeTimer.stop();

  int getWalkOnSpheresResultsSuccess = 
    getWalkOnSpheresResults(numWalksInProcess,			
			    parameters,
			    boundingSphere,
			    spheresNearestSurfacePointFinder,			
			    &threadRNGs,
			    &resultsZeno,
			    &walkTime,
			    &reduceTime);

  if (getWalkOnSpheresResultsSuccess != 0) {
    return getWalkOnSpheresResultsSuccess;
  }

  initializeTimer.start();

  SpheresInsideOutsideTester * spheresInsideOutsideTester = NULL;

  double surfacePreprocessTime = 0;

  initializeTimer.stop();

  int preprocessInteriorSuccess = 
    preprocessInterior(spheresLoaded,
		       parameters,
		       spheresModel,
		       &spheresInsideOutsideTester,
		       &surfacePreprocessTime);

  if (preprocessInteriorSuccess != 0) {
    return preprocessInteriorSuccess;
  }

  initializeTimer.start();

  ResultsInterior * resultsInterior = NULL;

  double sampleTime       = 0;
  double volumeReduceTime = 0;

  initializeTimer.stop();

  int getInteriorResultsSuccess = 
    getInteriorResults(numSamplesInProcess,			
		       parameters,
		       boundingSphere,
		       spheresInsideOutsideTester,			
		       &threadRNGs,
		       &resultsInterior,
		       &sampleTime,
		       &volumeReduceTime);

  if (getInteriorResultsSuccess != 0) {
    return getInteriorResultsSuccess;
  }

  initializeTime += initializeTimer.getTime();

  printOutput(boundingSphere,
	      resultsInterior,
	      resultsZeno,
	      parameters,
	      initializeTime,
	      readTime,
	      broadcastTime,
	      preprocessTime,
	      walkTime,
	      reduceTime,
	      surfacePreprocessTime,
	      sampleTime,
	      volumeReduceTime);

  savePointFiles(*resultsInterior,
		 *resultsZeno,
		 parameters);

  delete resultsZeno;
  delete resultsInterior;

  delete spheresNearestSurfacePointFinder;

  delete spheresInsideOutsideTester;

#ifdef USE_MPI
  MPI_Finalize();
#endif

  totalTimer.stop();

  if (mpiRank == 0) {
    std::cout << std::fixed
	      << "Total Time (s): " << totalTimer.getTime() << std::endl
	      << std::endl;

    printTime("End time: ");
  }

  return 0;
}

/// Get all input to the program, including both parameters and data, from the
/// command line and input files.
///
int
getInput(int argc, char **argv,
	 Parameters * parameters,
	 long long * numWalksInProcess,
	 long long * numSamplesInProcess,
	 Spheres<double> * spheres,
	 bool * spheresLoaded,
	 double * initializeTime,
	 double * readTime,
	 double * broadcastTime) {

  Timer initializeTimer;
  initializeTimer.start();

  parameters->parseCommandLine(argc, argv);

  if (parameters->getTotalNumWalksWasSet() &&
      (parameters->getMaxErrorCapacitanceWasSet() ||
       parameters->getMaxErrorPolarizabilityWasSet())) {

    std::cout << "Error: Cannot specify both number of walks and capacitance "
	      << "or polarizability tensor error" << std::endl;

    return 1;
  }

  if (parameters->getTotalNumSamplesWasSet() &&
      parameters->getMaxErrorVolumeWasSet()) {

    std::cout << "Error: Cannot specify both number of interior samples and "
	      << "volume error" << std::endl;

    return 1;
  }

  if (parameters->getComputeFormWasSet() &&
      !(parameters->getTotalNumSamplesWasSet() |
	parameters->getMaxErrorVolumeWasSet())) {

    std::cout << "Error: Must specify number of interior samples or volume "
	      << "error if computing form factors" << std::endl;

    return 1;
  }

  *numWalksInProcess = 
    computeNumInProcess(parameters->getMpiSize(), 
			parameters->getMpiRank(), 
			parameters->getTotalNumWalks());

  *numSamplesInProcess = 
    computeNumInProcess(parameters->getMpiSize(), 
			parameters->getMpiRank(), 
			parameters->getTotalNumSamples());

  if (parameters->getPrintBenchmarks() && 
      parameters->getMpiRank() == 0) {

    printRAM("RAM after initialization: ");
  }

  initializeTimer.stop();

  getBodData(parameters,
	     spheres,
	     readTime,
	     broadcastTime);

  initializeTimer.start();

  if (parameters->getPrintBenchmarks() && 
      parameters->getMpiRank() == 0) {

    printRAM("RAM after loading input data: ");
  }

  *spheresLoaded = !spheres->isEmpty();

  initializeTimer.stop();
  *initializeTime += initializeTimer.getTime();

  return 0;
}

/// Build the data structure used for the Walk-on-Spheres algorithm from either
/// voxel or sphere data.
///
int
preprocessWalkOnSpheres(bool spheresLoaded,
			Parameters const & parameters,
			Spheres<double> const & spheres,
			Sphere<double> * boundingSphere,
			SpheresModel * spheresModel,
			SpheresNearestSurfacePointFinder * *
			spheresNearestSurfacePointFinder,
			double * preprocessTime) {

  Timer preprocessTimer;
  preprocessTimer.start();

  if (spheresLoaded) {
    *boundingSphere = BoundingSphereGenerator::generate(spheres.getVector());

    spheresModel->preprocess(spheres.getVector(), 
			     parameters.getFracErrorBound());

    *spheresNearestSurfacePointFinder = 
      new SpheresNearestSurfacePointFinder(*spheresModel);
  }
  else {
    std::cout << "Warning: no spheres loaded" << std::endl;
  }

  if (parameters.getPrintBenchmarks() && 
      parameters.getMpiRank() == 0) {

    printRAM("RAM after building spatial data structure: ");
  }

  preprocessTimer.stop();
  *preprocessTime = preprocessTimer.getTime();

  return 0;
}

/// Allocate a random number generator for each thread, ensuring that each has
/// a unique stream ID across MPI processes.
///
void setupRNGs(Parameters const & parameters,
	       std::vector<RandomNumberGenerator> * threadRNGs) {

  int numStreams = 
    parameters.getNumThreads() * 
    parameters.getMpiSize();

  threadRNGs->reserve(parameters.getNumThreads());

  for (int threadNum = 0; threadNum < parameters.getNumThreads(); threadNum++) {
    int streamNum = 
      parameters.getMpiRank() * 
      parameters.getNumThreads() + 
      threadNum;

    threadRNGs->emplace_back(streamNum, numStreams, parameters.getSeed());
  }
}

/// Perform Walk-on-Spheres walks until the stopping condition is achieved
/// (either number of walks or error) and perform a parallel reduction on
/// the results.
///
int
getWalkOnSpheresResults(long long numWalksInProcess,			
			Parameters const & parameters,
			Sphere<double> const & boundingSphere,
			SpheresNearestSurfacePointFinder const * 
			spheresNearestSurfacePointFinder,			
			std::vector<RandomNumberGenerator> * threadRNGs,
			ResultsZeno * * resultsZeno,
			double * walkTime,
			double * reduceTime) {

  Timer reduceTimer;

  bool saveHitPoints = !parameters.getSurfacePointsFileName().empty();

  if (parameters.getTotalNumWalksWasSet()) {

    *resultsZeno = new ResultsZeno(boundingSphere,
				   parameters.getNumThreads(),
				   saveHitPoints);

    doWalkOnSpheresSelector(parameters,
			    numWalksInProcess,
			    boundingSphere, 
			    spheresNearestSurfacePointFinder,
			    threadRNGs,
			    *resultsZeno,
			    walkTime);

    reduceTimer.start();
    (*resultsZeno)->reduce();
    reduceTimer.stop();
  }
  else if (parameters.getMaxErrorCapacitanceWasSet() ||
	   parameters.getMaxErrorPolarizabilityWasSet()) {

    *resultsZeno = new ResultsZeno(boundingSphere,
				   parameters.getNumThreads(),
				   saveHitPoints);

    ResultsCompiler resultsCompiler(parameters);

    long long estimatedNumWalksRemaining = parameters.getMinTotalNumWalks();

    while (estimatedNumWalksRemaining > 0) {

      long long estimatedNumWalksRemainingInProcess = 
	computeNumInProcess(parameters.getMpiSize(), 
			    parameters.getMpiRank(), 
			    estimatedNumWalksRemaining);

      doWalkOnSpheresSelector(parameters,
			      estimatedNumWalksRemainingInProcess,
			      boundingSphere, 
			      spheresNearestSurfacePointFinder,
			      threadRNGs,
			      *resultsZeno,
			      walkTime);

      reduceTimer.start();
      (*resultsZeno)->reduce();
      reduceTimer.stop();

      resultsCompiler.compile(*resultsZeno,
			      NULL,
			      boundingSphere,
			      false);

      long long estimatedTotalNumWalks = 0;

      long long capacitanceEstimatedTotalNumWalks = 
  	estimateTotalNum(parameters.getMaxErrorCapacitance(),
  			 (*resultsZeno)->getNumWalks(),
  			 resultsCompiler.getCapacitance());

      estimatedTotalNumWalks = std::max(estimatedTotalNumWalks,
  					capacitanceEstimatedTotalNumWalks);

      long long polarizabilityEstimatedTotalNumWalks = 
  	estimateTotalNum(parameters.getMaxErrorPolarizability(),
  			 (*resultsZeno)->getNumWalks(),
  			 resultsCompiler.getMeanPolarizability());

      estimatedTotalNumWalks = std::max(estimatedTotalNumWalks,
  					polarizabilityEstimatedTotalNumWalks);

      estimatedNumWalksRemaining = 
	estimatedTotalNumWalks - (*resultsZeno)->getNumWalks();

      if (estimatedNumWalksRemaining > 0) {
	estimatedNumWalksRemaining = 
	  std::max(estimatedNumWalksRemaining,
		   parameters.getMinTotalNumWalks());
      }
    }
  }

  if (parameters.getPrintBenchmarks() && 
      parameters.getMpiRank() == 0) {

    printRAM("RAM after walk on spheres: ");
  }

  *reduceTime = reduceTimer.getTime();

  return 0;
}

/// Build the data structure used for the Interior Sampling algorithm from 
/// either voxel or sphere data.
///
int
preprocessInterior(bool spheresLoaded,
		   Parameters const & parameters,
		   SpheresModel const & spheresModel,
		   SpheresInsideOutsideTester * *
		   spheresInsideOutsideTester,
		   double * surfacePreprocessTime) {

  Timer surfacePreprocessTimer;
  surfacePreprocessTimer.start();

  if (spheresLoaded) {
    *spheresInsideOutsideTester = new SpheresInsideOutsideTester(spheresModel);
  }

  surfacePreprocessTimer.stop();
  *surfacePreprocessTime = surfacePreprocessTimer.getTime();

  return 0;
}

/// Perform Interior samples until the stopping condition is achieved
/// (either number of walks or error) and perform a parallel reduction on
/// the results.
///
int
getInteriorResults(long long numSamplesInProcess,			
		   Parameters const & parameters,
		   Sphere<double> const & boundingSphere,
		   SpheresInsideOutsideTester const * 
		   spheresInsideOutsideTester,			
		   std::vector<RandomNumberGenerator> * threadRNGs,
		   ResultsInterior * * resultsInterior,
		   double * sampleTime,
		   double * volumeReduceTime) {

  bool saveInteriorPoints = 
    !parameters.getInteriorPointsFileName().empty() ||
    parameters.getComputeFormWasSet();

  Timer volumeReduceTimer;

  if (parameters.getTotalNumSamplesWasSet()) {

    *resultsInterior = new ResultsInterior(parameters.getNumThreads(),
					   saveInteriorPoints);

    doInteriorSamplingSelector(parameters,
			       numSamplesInProcess,
			       boundingSphere,
			       spheresInsideOutsideTester,
			       threadRNGs,
			       *resultsInterior,
			       sampleTime);

    volumeReduceTimer.start();
    (*resultsInterior)->reduce();
    volumeReduceTimer.stop();
  }
  else if (parameters.getMaxErrorVolumeWasSet()) {

    *resultsInterior = new ResultsInterior(parameters.getNumThreads(),
					   saveInteriorPoints);

    ResultsCompiler resultsCompiler(parameters);

    long long estimatedNumSamplesRemaining = parameters.getMinTotalNumSamples();

    while (estimatedNumSamplesRemaining > 0) {

      long long estimatedNumSamplesRemainingInProcess = 
	computeNumInProcess(parameters.getMpiSize(), 
			    parameters.getMpiRank(), 
			    estimatedNumSamplesRemaining);

      doInteriorSamplingSelector(parameters,
				 estimatedNumSamplesRemainingInProcess,
				 boundingSphere, 
				 spheresInsideOutsideTester,
				 threadRNGs,
				 *resultsInterior,
				 sampleTime);

      volumeReduceTimer.start();
      (*resultsInterior)->reduce();
      volumeReduceTimer.stop();

      resultsCompiler.compile(NULL,
			      *resultsInterior,
			      boundingSphere,
			      false);

      long long estimatedTotalNumSamples = 
	estimateTotalNum(parameters.getMaxErrorVolume(),
			 (*resultsInterior)->getNumSamples(),
			 resultsCompiler.getVolume());

      estimatedNumSamplesRemaining = 
	estimatedTotalNumSamples - (*resultsInterior)->getNumSamples();

      if (estimatedNumSamplesRemaining > 0) {
	estimatedNumSamplesRemaining = 
	  std::max(estimatedNumSamplesRemaining,
		   parameters.getMinTotalNumSamples());
      }
    }
  }

  if (parameters.getComputeFormWasSet()) {
    volumeReduceTimer.start();
    (*resultsInterior)->gatherHitPoints();
    volumeReduceTimer.stop();
  }

  if (parameters.getPrintBenchmarks() && 
      parameters.getMpiRank() == 0) {

    printRAM("RAM after interior samples: ");
  }

  *volumeReduceTime = volumeReduceTimer.getTime();

  return 0;
}

/// Divides a total number of samples as evenly as possible between a set of
/// MPI processes of a certain size, and returns how many samples the MPI
/// process of the given rank is responsible for.
/// 
long long 
computeNumInProcess(int mpiSize, int mpiRank,
		    long long totalNumSamples) {

  long long numSamplesInProcess = totalNumSamples / mpiSize;

  if (mpiRank < totalNumSamples % mpiSize) {
    numSamplesInProcess ++;
  }

  return numSamplesInProcess;
}

/// Parses the bod file given as input and extracts sphere data, voxel data,
/// and parameters.
///
void
parseBodFile(Parameters * parameters, 
	     Spheres<double> * spheres) {

  std::string fileName = parameters->getInputFileName();

  std::ifstream inputFile;

  inputFile.open(fileName, std::ifstream::in);

  if (!inputFile.is_open()) {
    std::cout << "Error opening input file " << fileName << std::endl;
    exit(1);
  }

  Parser parser(inputFile, parameters, spheres);

  parser.parse();

  inputFile.close();
}

/// Gets the data from the bod file given as input either by parsing the file
/// or by an MPI brodcast, depending on the MPI rank of the process.
///
void
getBodData(Parameters * parameters,
	   Spheres<double> * spheres,
	   double * readTime,
	   double * broadcastTime) {

  Timer readTimer;
  Timer broadcastTimer;

  if (parameters->getMpiRank() == 0) {
    readTimer.start();
    parseBodFile(parameters, spheres);
    readTimer.stop();

    broadcastTimer.start();
    spheres->mpiSend();
    parameters->mpiSend();
    broadcastTimer.stop();
  }
  else {
    broadcastTimer.start();
    spheres->mpiReceive();
    parameters->mpiReceive();
    broadcastTimer.stop();
  }

  *readTime      = readTimer.getTime();
  *broadcastTime = broadcastTimer.getTime();
}

/// Sets the default values for parameters that have them if the parameters
/// have not already been set.
///
void
computeDefaultParameters(Parameters * parameters,
			 Sphere<double> * boundingSphere) {

  const double defaultSkinThicknessFactor = 0.000001;

  if (parameters->getLaunchCenterWasSet()) {
    boundingSphere->setCenter(parameters->getLaunchCenter());
  }
  else {
    parameters->setLaunchCenter(boundingSphere->getCenter());
  }

  if (parameters->getLaunchRadiusWasSet()) {
    boundingSphere->setRadius(parameters->getLaunchRadius());
  }
  else {
    parameters->setLaunchRadius(boundingSphere->getRadius());
  }

  if (!parameters->getSkinThicknessWasSet()) {
    parameters->setSkinThickness(boundingSphere->getRadius() * 
				 defaultSkinThicknessFactor);
  }

  if (!parameters->getLengthScaleWasSet()) {
    parameters->setLengthScale(1, Units::Length::L);
  }
}

/// Estimates the total number of samples that will be required to acheive the
/// requested error, given the number of samples taken so far and the current
/// value and error.
///
long long 
estimateTotalNum(double requestedError,
		 long long numSoFar,
		 Uncertain<double> const & currentValue) {

  double requestedVariance = 
    pow((requestedError / 100) * currentValue.getMean(), 2);

  long long estimatedTotalNum = 
    ceil(currentValue.getVariance() * numSoFar / requestedVariance);

  return estimatedTotalNum;
}

/// Makes a call to the doWalkOnSpheres function with either sphere or voxel
/// data, depending which is not NULL.
///
void 
doWalkOnSpheresSelector(Parameters const & parameters,
			long long numWalksInProcess,
			Sphere<double> const & boundingSphere, 
			SpheresNearestSurfacePointFinder const * 
			spheresNearestSurfacePointFinder,
			std::vector<RandomNumberGenerator> * threadRNGs,
			ResultsZeno * resultsZeno,
			double * walkTime) {

  if (spheresNearestSurfacePointFinder != NULL) {
    doWalkOnSpheres(parameters.getNumThreads(),
		    numWalksInProcess,
		    boundingSphere, 
		    *spheresNearestSurfacePointFinder,
		    parameters.getFracErrorBound(),
		    parameters.getSkinThickness(),
		    threadRNGs,
		    resultsZeno,
		    walkTime);
  }
  else {
    assert(0);
  }
}

/// Launches a set of Walk-on-Spheres walks in each of a set of parallel 
/// threads.
///
template <class NearestSurfacePointFinder>
void
doWalkOnSpheres(int numThreads,
		long long numWalksInProcess,
		Sphere<double> const & boundingSphere, 
		NearestSurfacePointFinder const & 
		nearestSurfacePointFinder,
		double fracErrorBound,
		double shellThickness,
		std::vector<RandomNumberGenerator> * threadRNGs,
		ResultsZeno * resultsZeno,
		double * walkTime) {

  Timer walkTimer;
  walkTimer.start();

  std::thread * * threads = new std::thread *[numThreads];

  for (int threadNum = 0; threadNum < numThreads; threadNum++) {

    long long numWalksInThread = numWalksInProcess / numThreads;

    if (threadNum < numWalksInProcess % numThreads) {
      numWalksInThread ++;
    }

    threads[threadNum] = 
      new std::thread(doWalkOnSpheresThread<NearestSurfacePointFinder>,
		      boundingSphere, 
		      nearestSurfacePointFinder,
		      threadNum,
		      fracErrorBound,
		      shellThickness,
		      numWalksInThread,
		      &(threadRNGs->at(threadNum)),
		      resultsZeno);
  }

  for (int threadNum = 0; threadNum < numThreads; threadNum++) {
    threads[threadNum]->join();
  }

  for (int threadNum = 0; threadNum < numThreads; threadNum++) {
    delete threads[threadNum];
  }

  delete [] threads;

  walkTimer.stop();
  *walkTime += walkTimer.getTime();
}

/// Launches a given number of Walk-on-Spheres walks and records the results.
/// Runs in a single thread.
///
template <class NearestSurfacePointFinder>
void
doWalkOnSpheresThread(Sphere<double> const & boundingSphere, 
		      NearestSurfacePointFinder const & 
		      nearestSurfacePointFinder,
		      int threadNum,
		      double fracErrorBound,
		      double shellThickness,
		      long long numWalks,
		      RandomNumberGenerator * randomNumberGenerator,
		      ResultsZeno * resultsZeno) {

  WalkerExterior<double, 
		 RandomNumberGenerator,
		 NearestSurfacePointFinder,
		 RandomSpherePointGenerator,
		 BiasedSpherePointGenerator>
    walker(randomNumberGenerator, 
	   boundingSphere, 
           nearestSurfacePointFinder,
           fracErrorBound,
	   shellThickness);

  for (long long walkNum = 0; walkNum < numWalks; walkNum++) {

    bool hitObject = false;
    int numSteps   = 0;

    Vector3<double> startPoint;
    Vector3<double> endPoint;
    Vector3<double> normal;

    walker.walk(&hitObject, &numSteps,
		&startPoint, &endPoint, &normal);

    if (hitObject) {
      resultsZeno->recordHit(threadNum, 
			     startPoint, endPoint, normal,
			     randomNumberGenerator);
    }
    else {
      resultsZeno->recordMiss(threadNum);
    }
  }
}

/// Makes a call to the doInteriorSampling function with either sphere or voxel
/// data, depending which is not NULL.
///
void
doInteriorSamplingSelector(Parameters const & parameters,
			   long long numSamplesInProcess,
			   Sphere<double> const & boundingSphere,
			   SpheresInsideOutsideTester const *  
			   spheresInsideOutsideTester,
			   std::vector<RandomNumberGenerator> * threadRNGs,
			   ResultsInterior * resultsInterior,
			   double * sampleTime) {

  if (spheresInsideOutsideTester != NULL) {
    doInteriorSampling(parameters.getNumThreads(),
		       numSamplesInProcess,
		       boundingSphere, 
		       *spheresInsideOutsideTester,
		       parameters.getFracErrorBound(),
		       threadRNGs,
		       resultsInterior,
		       sampleTime);
  }
  else {
    assert(0);
  }
}

/// Launches a set of Interior samples in each of a set of parallel 
/// threads.
///
template <class InsideOutsideTester>
void
doInteriorSampling(int numThreads,
		   long long numSamplesInProcess,
		   Sphere<double> const & boundingSphere, 
		   InsideOutsideTester const & insideOutsideTester,
		   double fracErrorBound,
		   std::vector<RandomNumberGenerator> * threadRNGs,
		   ResultsInterior * resultsInterior,
		   double * sampleTime) {

  Timer sampleTimer;
  sampleTimer.start();

  std::thread * * threads = new std::thread *[numThreads];

  for (int threadNum = 0; threadNum < numThreads; threadNum++) {

    long long numSamplesInThread = numSamplesInProcess / numThreads;

    if (threadNum < numSamplesInProcess % numThreads) {
      numSamplesInThread ++;
    }

    threads[threadNum] = 
      new std::thread(doInteriorSamplingThread<InsideOutsideTester>,
		      boundingSphere, 
		      insideOutsideTester,
		      threadNum,
		      fracErrorBound,
		      numSamplesInThread,
		      &(threadRNGs->at(threadNum)),
		      resultsInterior);
  }

  for (int threadNum = 0; threadNum < numThreads; threadNum++) {
    threads[threadNum]->join();
  }

  for (int threadNum = 0; threadNum < numThreads; threadNum++) {
    delete threads[threadNum];
  }

  delete [] threads;

  sampleTimer.stop();
  *sampleTime += sampleTimer.getTime();
}

/// Performs a given number of Interior samples and records the results.
/// Runs in a single thread.
///
template <class InsideOutsideTester>
void
doInteriorSamplingThread(Sphere<double> const & boundingSphere, 
			 InsideOutsideTester const & insideOutsideTester,
			 int threadNum,
			 double fracErrorBound,
			 long long numSamples,
			 RandomNumberGenerator * randomNumberGenerator,
			 ResultsInterior * resultsInterior) {

  SamplerInterior<double, 
		 RandomNumberGenerator,
		 InsideOutsideTester,
		 RandomBallPointGenerator>
    sampler(randomNumberGenerator, 
	    boundingSphere, 
	    insideOutsideTester,
	    fracErrorBound);

  for (long long sampleNum = 0; sampleNum < numSamples; sampleNum++) {

    bool hitObject = false;

    Vector3<double> hitPoint;

    sampler.sample(&hitObject,
		   &hitPoint);

    if (hitObject) {
      resultsInterior->recordHit(threadNum,
				 hitPoint);
    }
    else {
      resultsInterior->recordMiss(threadNum);
    }
  }
}

/// Prints parameters, results, and (optionally) detailed running time
/// benchmarks on MPI process 0.
///
void
printOutput(Sphere<double> const & boundingSphere,
	    ResultsInterior const * resultsInterior,
	    ResultsZeno const * resultsZeno,
	    Parameters const & parameters, 
	    double initializeTime,
	    double readTime,
	    double broadcastTime,
	    double preprocessTime,
	    double walkTime,
	    double reduceTime,
	    double surfacePreprocessTime,
	    double sampleTime,
	    double volumeReduceTime) {

  if (parameters.getMpiRank() == 0) {
    std::cout << std::endl
	      << "Parameters" << std::endl
	      << "----------" << std::endl
	      << std::endl;

    parameters.print();

    std::cout << std::endl
	      << "Results" << std::endl
	      << "-------" << std::endl
	      << std::endl;

    if (resultsZeno != NULL) {

      std::cout << "Number of walks performed: " 
		<< (long long)resultsZeno->getNumWalks() << std::endl
		<< std::endl;
    }

    if (resultsInterior != NULL) {

      std::cout << "Number of interior samples taken: " 
		<< (long long)resultsInterior->getNumSamples() << std::endl
		<< std::endl;
    }

    ResultsCompiler resultsCompiler(parameters);

    resultsCompiler.compile(resultsZeno,
			    resultsInterior,
			    boundingSphere,
			    parameters.getComputeFormWasSet());

    resultsCompiler.print(parameters.getPrintCounts());

    if (parameters.getPrintBenchmarks()) {
      std::cout << std::fixed
		<< "Initialize (s):         " << initializeTime << std::endl
		<< "Read (s):               " << readTime << std::endl
		<< "Broadcast (s):          " << broadcastTime << std::endl
		<< "Centers Preprocess (s): " << preprocessTime << std::endl
		<< "Exterior Walk (s):      " << walkTime << std::endl
		<< "Exterior Reduce (s):    " << reduceTime << std::endl
		<< "Surface Preprocess (s): " << surfacePreprocessTime 
		<< std::endl
		<< "Volume Sample (s):      " << sampleTime << std::endl
		<< "Volume Reduce (s):      " << volumeReduceTime << std::endl
		<< std::endl;
    }
  }
}

/// Writes Walk-on-Spheres and Interior Sampling hit points to disk.
/// 
void
savePointFiles(ResultsInterior & resultsInterior,
	       ResultsZeno & resultsZeno,
	       Parameters const & parameters) {

  if (!parameters.getSurfacePointsFileName().empty()) {
    resultsZeno.gatherHitPoints();

    writePoints(parameters.getSurfacePointsFileName(), 
		resultsZeno.getPoints(), 
		resultsZeno.getCharges());
  }

  if (!parameters.getInteriorPointsFileName().empty()) {
    resultsInterior.gatherHitPoints();

    writePoints(parameters.getInteriorPointsFileName(), 
		resultsInterior.getPoints(), 
		NULL);
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
printRAM(std::string const & label) {
  std::ifstream statusFile("/proc/self/status");

  std::string line;

  while (statusFile.good()) {
    std::getline(statusFile, line);

    size_t pos = line.find("VmRSS:");

    if (pos != std::string::npos) {
      std::cout << label << line.substr(pos + 6) << std::endl;
    }
  }

  statusFile.close();
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

    std::cout << "Error opening output file " << fileName << std::endl;
    exit(1);
  }

  for (unsigned int i = 0; i < points->size(); i++) {

    if (charges != NULL) {
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

// ================================================================

// Local Variables:
// time-stamp-line-limit: 30
// End:
