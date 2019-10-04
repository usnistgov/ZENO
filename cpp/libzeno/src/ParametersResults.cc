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

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "ParametersResults.h"

// ================================================================

using namespace zeno;

ParametersResults::ParametersResults() 
  : computeForm(false),
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

ParametersResults::~ParametersResults() {

}

void
ParametersResults::setComputeForm(bool computeForm) {
  this->computeForm = computeForm;
}

bool
ParametersResults::getComputeForm() const {
  return computeForm;
}

void
ParametersResults::setLengthScale(double number, Units::Length unit) {
  lengthScale     = number;
  lengthScaleUnit = unit;

  lengthScaleWasSet = true;
}

double
ParametersResults::getLengthScaleNumber() const {
  return lengthScale;
}

Units::Length
ParametersResults::getLengthScaleUnit() const {
  return lengthScaleUnit;
}

bool
ParametersResults::getLengthScaleWasSet() const {
  return lengthScaleWasSet;
}

void 
ParametersResults::setTemperature(double number, Units::Temperature unit) {
  temperature     = number;
  temperatureUnit = unit;

  temperatureWasSet = true;
}

double
ParametersResults::getTemperatureNumber() const {
  return temperature;
}

Units::Temperature
ParametersResults::getTemperatureUnit() const {
  return temperatureUnit;
}

bool
ParametersResults::getTemperatureWasSet() const {
  return temperatureWasSet;
}

void
ParametersResults::setMass(double number, Units::Mass unit) {
  mass     = number;
  massUnit = unit;

  massWasSet = true;
}

double 
ParametersResults::getMassNumber() const {
  return mass;
}

Units::Mass
ParametersResults::getMassUnit() const {
  return massUnit;
}

bool
ParametersResults::getMassWasSet() const {
  return massWasSet;
}

void 
ParametersResults::setSolventViscosity(double number, Units::Viscosity unit) {
  solventViscosity     = number;
  solventViscosityUnit = unit;

  solventViscosityWasSet = true;
}

double
ParametersResults::getSolventViscosityNumber() const {
  return solventViscosity;
}

Units::Viscosity
ParametersResults::getSolventViscosityUnit() const {
  return solventViscosityUnit;
}

bool
ParametersResults::getSolventViscosityWasSet() const {
  return solventViscosityWasSet;
}

void 
ParametersResults::setBuoyancyFactor(double buoyancyFactor) {
  this->buoyancyFactor = buoyancyFactor;

  buoyancyFactorWasSet = true;
}

double
ParametersResults::getBuoyancyFactor() const {
  return buoyancyFactor;
}

bool
ParametersResults::getBuoyancyFactorWasSet() const {
  return buoyancyFactorWasSet;
}

void 
ParametersResults::mpiBroadcast(int root) {
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

/// Broadcasts the parameters that can be set in the bod file over MPI.
///
void 
ParametersResults::serializeMpiBroadcast(int root) const {
#ifdef USE_MPI
  const int countToSend = 15;
  
  double parametersArray[countToSend];
  
  parametersArray[0] = getLengthScaleNumber();
  parametersArray[1] = (double)getLengthScaleUnit();
  parametersArray[2] = (double)getLengthScaleWasSet();
  parametersArray[3] = getTemperatureNumber();
  parametersArray[4] = (double)getTemperatureUnit();
  parametersArray[5] = (double)getTemperatureWasSet();
  parametersArray[6] = getMassNumber();
  parametersArray[7] = (double)getMassUnit();
  parametersArray[8] = (double)getMassWasSet();
  parametersArray[9] = getSolventViscosityNumber();
  parametersArray[10] = (double)getSolventViscosityUnit();
  parametersArray[11] = (double)getSolventViscosityWasSet();
  parametersArray[12] = getBuoyancyFactor();
  parametersArray[13] = (double)getBuoyancyFactorWasSet();
  parametersArray[14] = (double)getComputeForm();

  MPI_Bcast(parametersArray, countToSend, MPI_DOUBLE, root, MPI_COMM_WORLD);
#endif
}

/// Receives the parameters that can be set in the bod file over MPI.
///
void 
ParametersResults::mpiBroadcastDeserialize(int root) {
#ifdef USE_MPI
  const int countToReceive = 15;
  
  double parametersArray[countToReceive];

  MPI_Bcast(parametersArray, countToReceive, MPI_DOUBLE, root, MPI_COMM_WORLD);

  if ((bool)parametersArray[2]) {
    setLengthScale(parametersArray[0], 
		   (Units::Length)parametersArray[1]);
  }

  if ((bool)parametersArray[5]) {
    setTemperature(parametersArray[3],
		   (Units::Temperature)parametersArray[4]);
  }

  if ((bool)parametersArray[8]) {
    setMass(parametersArray[6],
	    (Units::Mass)parametersArray[7]);
  }

  if ((bool)parametersArray[11]) {
    setSolventViscosity(parametersArray[9],
			(Units::Viscosity)parametersArray[10]);
  }

  if ((bool)parametersArray[13]) {
    setBuoyancyFactor(parametersArray[12]);
  }

  setComputeForm((bool)parametersArray[14]);
#endif
}
