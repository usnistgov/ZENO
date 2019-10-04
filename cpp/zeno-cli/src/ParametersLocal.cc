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

#include "ParametersLocal.h"

// ================================================================

ParametersLocal::ParametersLocal() 
  : inputFileName(),
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
