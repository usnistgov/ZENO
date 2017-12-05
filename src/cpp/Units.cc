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
// Author:  Derek Juba <derek.juba@nist.gov>
// Date:    Tue Jan 12 17:37:19 2016 EDT
// 
// Time-stamp: <2017-12-05 16:37:03 dcj>
// 
// ================================================================

#include <cassert>
#include <cmath>

#include "Units.h"

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
    return pow(10, -7);
  case Length::A:
    return pow(10, -8);
  default:
    assert(0);
  }
}
 
Uncertain<double> 
Units::getFactor(Mass fromUnit, Mass toUnit) {
  assert(toUnit == Mass::g);

  switch(fromUnit) {
  case Mass::Da:
    return Uncertain<double>(6.02214*pow(10, 23), pow(10, 18 * 2));
  case Mass::kDa:
    return Uncertain<double>(6.02214*pow(10, 20), pow(10, 15 * 2));
  case Mass::g:
    return 1;
  case Mass::kg:
    return pow(10, -3);
  default:
    assert(0);
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
  } 
}

/// Returns Boltzmann's constant.
///
Uncertain<double> 
Units::kB() {
  return Uncertain<double>(1.38065*pow(10, -23), pow(10, -28 * 2));
}

// ================================================================

// Local Variables:
// time-stamp-line-limit: 30
// End:
