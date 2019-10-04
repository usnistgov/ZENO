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
// Created: Tue Jan 12 17:37:19 2016 EDT
// 
// ================================================================

#include <cassert>
#include <cmath>

#include "Units.h"

using namespace zeno;

// ================================================================

Uncertain<double> 
Units::getFactor(Length fromUnit, Length toUnit) {
  assert(toUnit == Length::cm);

  switch(fromUnit) {
  case Length::m: 
    return 100; 
  case Length::cm:
    return 1;
  case Length::nm:
    return std::pow(10, -7);
  case Length::A:
    return std::pow(10, -8);
  default:
    assert(0);
    exit(EXIT_FAILURE);
  }
}
 
Uncertain<double> 
Units::getFactor(Mass fromUnit, Mass toUnit) {
  assert(toUnit == Mass::g);

  switch(fromUnit) {
  case Mass::Da:
    return Uncertain<double>(6.02214*std::pow(10, 23), std::pow(10, 18 * 2));
  case Mass::kDa:
    return Uncertain<double>(6.02214*std::pow(10, 20), std::pow(10, 15 * 2));
  case Mass::g:
    return 1;
  case Mass::kg:
    return std::pow(10, -3);
  default:
    assert(0);
    exit(EXIT_FAILURE);
  }
}

Uncertain<double> 
Units::getFactor(Viscosity fromUnit, Viscosity toUnit) {
  assert(toUnit == Viscosity::cp);

  switch(fromUnit) {
  case Viscosity::cp:
    return 1;
  case Viscosity::p:
    return 100;
  default:
    assert(0);
    exit(EXIT_FAILURE);
  }
}

Uncertain<double> 
Units::getOffset(Temperature fromUnit, Temperature toUnit) {
  assert(toUnit == Temperature::K);

  switch(fromUnit) {
  case Temperature::C:
    return 273.15;
  case Temperature::K:
    return 0;
  default:
    assert(0);
    exit(EXIT_FAILURE);
  }
}

std::string 
Units::getName(Length unit) {
  switch(unit) {
  case Length::m:
    return "m";
  case Length::cm:
    return "cm";
  case Length::nm:
    return "nm";
  case Length::A:
    return "A";
  case Length::L:
    return "L";
  default:
    assert(0);
    return "*UNKNOWN*";
  }
}

std::string 
Units::getName(Temperature unit) {
  switch(unit) {
  case Temperature::C:
    return "C";
  case Temperature::K:
    return "K";
  default:
    assert(0);
    return "*UNKNOWN*";
  }
}

std::string 
Units::getName(Mass unit) {
  switch(unit) {
  case Mass::Da:
    return "Da";
  case Mass::kDa:
    return "kDa";
  case Mass::g:
    return "g";
  case Mass::kg:
    return "kg";
  default:
    assert(0);
    return "*UNKNOWN*";
  }
}

std::string 
Units::getName(Viscosity unit) {
  switch(unit) {
  case Viscosity::p:
    return "p";
  case Viscosity::cp:
    return "cp";
  default:
    assert(0);
    return "*UNKNOWN*";
  } 
}

/// Returns Boltzmann's constant.
///
Uncertain<double> 
Units::kB() {
  return Uncertain<double>(1.38065*std::pow(10, -23), std::pow(10, -28 * 2));
}

