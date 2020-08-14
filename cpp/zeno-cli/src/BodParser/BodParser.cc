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
// Created: Mon Feb 24 11:24:59 2014 EDT
//
// ================================================================

#include "BodParser.h"

#include "Geometry/Vector3.h"

// ================================================================

using namespace zeno;

bod_parser::BodParser::BodParser
(ParametersLocal * parametersLocal,
 std::istream &in,
 ParametersWalkOnSpheres * parametersWalkOnSpheres,
 ParametersInteriorSampling * parametersInteriorSampling,
 ParametersResults * parametersResults,
 MixedModel<double> * model) :
  d_scanner(in),
  parametersWalkOnSpheres(parametersWalkOnSpheres),
  parametersInteriorSampling(parametersInteriorSampling),
  parametersResults(parametersResults),
  parametersLocal(parametersLocal),
  model(model) {

}

void bod_parser::BodParser::addSphere(double x, double y, double z, double r) {
  if (r < 0) {
    std::cerr << "Error: tried to add sphere with negative radius " << r
	      << std::endl;

    exit(1);
  }
  
  Vector3<double> center(x, y, z);

  Sphere<double> sphere(center, r);

  model->addSphere(sphere);
}

void bod_parser::BodParser::addCube(double x, double y, double z, double s) {
  if (s < 0) {
    std::cerr << "Error: tried to add cube with negative offset " << s
	      << std::endl;

    exit(1);
  }
  
  Vector3<double> minCoords(x,     y,     z);
  Vector3<double> maxCoords(x + s, y + s, z + s);

  Cuboid<double> cuboid(minCoords, maxCoords);
  
  model->addCuboid(cuboid);
}

void bod_parser::BodParser::addCuboid(double x1, double y1, double z1,
				      double x2, double y2, double z2) {
  
  if (x2 < x1) std::swap(x1, x2);
  if (y2 < y1) std::swap(y1, y2);
  if (z2 < z1) std::swap(z1, z2);
  
  Vector3<double> minCoords(x1, y1, z1);
  Vector3<double> maxCoords(x2, y2, z2);

  Cuboid<double> cuboid(minCoords, maxCoords);
  
  model->addCuboid(cuboid);
}

void bod_parser::BodParser::addVoxels(std::string voxelsFileName) {
  const std::string inputFilePath = parametersLocal->getInputFileName();

  size_t directoryEndPos = inputFilePath.find_last_of('/');

  if (directoryEndPos == std::string::npos) {
    directoryEndPos = 0;
  }
  else {
    directoryEndPos += 1;
  }
  
  std::string voxelsFilePath = inputFilePath;
  
  voxelsFilePath.replace(directoryEndPos, std::string::npos, voxelsFileName);
  
  bool voxelsLoaded = 
    model->loadVoxels(voxelsFilePath.c_str());

  if (!voxelsLoaded) {
    std::cerr << "Error loading voxel file: " << voxelsFileName 
	      << std::endl;
      
    exit(1);
  }
}

void bod_parser::BodParser::addTrajectory(std::string xyzFileName, std::string mapFileName) {
  parametersLocal->setXyzInputFileName(xyzFileName);

  parametersLocal->setMapInputFileName(mapFileName);
}

void bod_parser::BodParser::setST(double skinThickness) {
  parametersWalkOnSpheres->setSkinThickness(skinThickness);
}

void bod_parser::BodParser::setRLAUNCH(double launchRadius) {
  parametersWalkOnSpheres->setLaunchRadius(launchRadius);

  parametersInteriorSampling->setLaunchRadius(launchRadius);
}

void bod_parser::BodParser::setHUNITS(double number, std::string unitString) {
  Units::Length units;

  if (unitString == "m") {
    units = Units::Length::m;
  }
  else if (unitString == "cm") {
    units = Units::Length::cm;
  }
  else if (unitString == "nm") {
    units = Units::Length::nm;
  }
  else if (unitString == "A") {
    units = Units::Length::A;
  }
  else if (unitString == "L") {
    units = Units::Length::L;
  }
  else {
    std::cerr << "Error: invalid length unit " << unitString << std::endl;
    exit(1);
  }

  parametersResults->setLengthScale(number, units);
}

void bod_parser::BodParser::setUNITS(std::string unitString) {
  setHUNITS(1, unitString);
}

void bod_parser::BodParser::setTEMP(double number, std::string unitString) {
  Units::Temperature units;

  if (unitString == "C") {
    units = Units::Temperature::C;
  }
  else if (unitString == "K") {
    units = Units::Temperature::K;
  }
  else {
    std::cerr << "Error: invalid temperature unit " << unitString << std::endl;
    exit(1);
  }

  parametersResults->setTemperature(number, units);
}

void bod_parser::BodParser::setMASS(double number, std::string unitString) {
  Units::Mass units;

  if (unitString == "Da") {
    units = Units::Mass::Da;
  }
  else if (unitString == "kDa") {
    units = Units::Mass::kDa;
  }
  else if (unitString == "g") {
    units = Units::Mass::g;
  }
  else if (unitString == "kg") {
    units = Units::Mass::kg;
  }
  else {
    std::cerr << "Error: invalid mass unit " << unitString << std::endl;
    exit(1);
  }

  parametersResults->setMass(number, units);
}

void bod_parser::BodParser::setVISCOSITY(double number,
					 std::string unitString) {
  Units::Viscosity units;

  if (unitString == "p") {
    units = Units::Viscosity::p;
  }
  else if (unitString == "cp") {
    units = Units::Viscosity::cp;
  }
  else {
    std::cerr << "Error: invalid viscosity unit " << unitString << std::endl;
    exit(1);
  }

  parametersResults->setSolventViscosity(number, units);
}

void bod_parser::BodParser::setBF(double buoyancyFactor) {
  parametersResults->setBuoyancyFactor(buoyancyFactor);
}

