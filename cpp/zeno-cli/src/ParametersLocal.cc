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
// Created: 2019-05-30
// 
// ================================================================

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "ParametersLocal.h"

// ================================================================

ParametersLocal::ParametersLocal() 
  : inputFileName(),
    xyzInputFileName(),
    xyzInputFileNameWasSet(false),
    mapInputFileName(),
    mapInputFileNameWasSet(false),
    csvOutputFileName(),
    csvOutputFileNameWasSet(false),
    mpiSize(1),
    mpiRank(0),
    surfacePointsFileName(""),
    interiorPointsFileName(""),
    printCounts(),
    printBenchmarks() {

}

ParametersLocal::~ParametersLocal() {

}

void
ParametersLocal::setInputFileName(std::string const & inputFileName) {
  this->inputFileName = inputFileName;
}

std::string 
ParametersLocal::getInputFileName() const { 
  return inputFileName;
}

void
ParametersLocal::setXyzInputFileName(std::string const & xyzInputFileName) {
  this->xyzInputFileName = xyzInputFileName;

  xyzInputFileNameWasSet = true;
}

std::string
ParametersLocal::getXyzInputFileName() const {
  return xyzInputFileName;
}

bool
ParametersLocal::getXyzInputFileNameWasSet() const {
  return xyzInputFileNameWasSet;
}

void
ParametersLocal::setMapInputFileName(std::string const & mapInputFileName) {
  this->mapInputFileName = mapInputFileName;

  mapInputFileNameWasSet = true;
}

std::string
ParametersLocal::getMapInputFileName() const {
  return mapInputFileName;
}

bool
ParametersLocal::getMapInputFileNameWasSet() const {
  return mapInputFileNameWasSet;
}

void
ParametersLocal::setCsvOutputFileName(std::string const & csvOutputFileName) {
  this->csvOutputFileName = csvOutputFileName;

  csvOutputFileNameWasSet = true;
}

std::string
ParametersLocal::getCsvOutputFileName() const {
  return csvOutputFileName;
}

bool
ParametersLocal::getCsvOutputFileNameWasSet() const {
  return csvOutputFileNameWasSet;
}

int
ParametersLocal::getMpiSize() const {
  return mpiSize;
}

void
ParametersLocal::setMpiSize(int mpiSize) {
  this->mpiSize = mpiSize;
}

int
ParametersLocal::getMpiRank() const {
  return mpiRank;
}

void
ParametersLocal::setMpiRank(int mpiRank) {
  this->mpiRank = mpiRank;
}

void
ParametersLocal::setSurfacePointsFileName
(std::string const & surfacePointsFileName) {
  this->surfacePointsFileName = surfacePointsFileName;
}

std::string 
ParametersLocal::getSurfacePointsFileName() const {
  return surfacePointsFileName;
}

void
ParametersLocal::setInteriorPointsFileName
(std::string const & interiorPointsFileName) {
  this->interiorPointsFileName = interiorPointsFileName;
}

std::string 
ParametersLocal::getInteriorPointsFileName() const {
  return interiorPointsFileName;
}

void
ParametersLocal::setPrintCounts(bool printCounts) {
  this->printCounts = printCounts;
}

bool 
ParametersLocal::getPrintCounts() const {
  return printCounts;
}

void
ParametersLocal::setPrintBenchmarks(bool printBenchmarks) {
  this->printBenchmarks = printBenchmarks;
}

bool 
ParametersLocal::getPrintBenchmarks() const {
  return printBenchmarks;
}

void                                                                                                                   ParametersLocal::mpiBroadcast(int root) {
  #ifdef USE_MPI
  int mpiSize = 1, mpiRank = 0;

  MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

  if (mpiSize > 1) {
    if (mpiRank == root) {
      serializeMpiBroadcast(root);
    }
    else {
      mpiBroadcastDeserialize(root);
    }
  }
  #endif
}

/// Broadcasts a subset of parameters over MPI.
///
void
ParametersLocal::serializeMpiBroadcast(int root) const {
#ifdef USE_MPI
  const int numLongLongsToSend = 2;

  long long longLongsArray[numLongLongsToSend];
  longLongsArray[0] = (long long)getXyzInputFileNameWasSet();
  longLongsArray[1] = (long long)getMapInputFileNameWasSet();
 
  MPI_Bcast(longLongsArray, numLongLongsToSend, MPI_LONG_LONG,
	    root, MPI_COMM_WORLD);
#endif
}

/// Receives a subset of parameters over MPI.
///
void
ParametersLocal::mpiBroadcastDeserialize(int root) {
#ifdef USE_MPI
  const int numLongLongsToReceive = 2;

  long long longLongsArray[numLongLongsToReceive];
  MPI_Bcast(longLongsArray, numLongLongsToReceive, MPI_LONG_LONG,
	    root, MPI_COMM_WORLD);

  if ((bool)longLongsArray[0]) {
    xyzInputFileNameWasSet = true;
  }

  if ((bool)longLongsArray[1]) {
    mapInputFileNameWasSet = true;
  }
#endif
}
