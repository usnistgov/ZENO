// ================================================================
/// 
/// Disclaimer:  IMPORTANT:  This software was developed at the
/// National Institute of Standards and Technology by employees of the
/// Federal Government in the course of their official duties.
/// Pursuant to title 17 Section 105 of the United States Code this
/// software is not subject to copyright protection and is in the
/// public domain.  This is an experimental system.  NIST assumes no
/// responsibility whatsoever for its use by other parties, and makes
/// no guarantees, expressed or implied, about its quality,
/// reliability, or any other characteristic.  We would appreciate
/// acknowledgement if the software is used.  This software can be
/// redistributed and/or modified freely provided that any derivative
/// works bear some notice that they are derived from it, and any
/// modified versions bear some notice that they have been modified.
/// 
/// ================================================================

// ================================================================
// 
// Authors: Derek Juba <derek.juba@nist.gov>
// Date:    Mon Feb 24 11:24:59 2014 EDT
//
// Time-stamp: <2016-09-22 13:10:45 dcj>
//
// ================================================================

#include "Parser.h"

#include "Geometry/Vector3.h"
#include "Geometry/Sphere.h"

// ================================================================

Parser::Parser(std::istream &in,
	       Parameters * parameters,
	       Spheres<double> * spheres) :
  d_scanner(in),
  parameters(parameters),
  spheres(spheres){

}

void Parser::addSphere(double x, double y, double z, double r) {
  Vector3<double> center(x, y, z);

  Sphere<double> sphere(center, r);

  spheres->add(sphere);
}

void Parser::setST(double skinThickness) {
  parameters->setSkinThickness(skinThickness);
}

void Parser::setRLAUNCH(double launchRadius) {
  parameters->setLaunchRadius(launchRadius);
}

void Parser::setHUNITS(double number, std::string unitString) {
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

  parameters->setLengthScale(number, units);
}

void Parser::setUNITS(std::string unitString) {
  setHUNITS(1, unitString);
}

void Parser::setTEMP(double number, std::string unitString) {
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

  parameters->setTemperature(number, units);
}

void Parser::setMASS(double number, std::string unitString) {
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

  parameters->setMass(number, units);
}

void Parser::setVISCOSITY(double number, std::string unitString) {
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

  parameters->setSolventViscosity(number, units);
}

void Parser::setBF(double buoyancyFactor) {
  parameters->setBuoyancyFactor(buoyancyFactor);
}

// ================================================================

// Local Variables:
// time-stamp-line-limit: 30
// End:
