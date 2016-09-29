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
// Date:    Tue Jan 12 17:05:37 2016 EDT
//
// Time-stamp: <2016-08-30 16:13:19 dcj>
//
// ================================================================

#ifndef UNITS_H_
#define UNITS_H_

// ================================================================

#include <string>

#include "Uncertain.h"

// ================================================================

/// Contains symbolic constants representing units.  Also contains functions 
/// for getting conversion factors and offsets and human-readable names.
///
class Units
{
public:
  enum class Length {m, cm, nm, A, L};
  enum class Temperature {C, K};
  enum class Mass {Da, kDa, g, kg};
  enum class Viscosity {p, cp};

  static Uncertain<double> getFactor(Length fromUnit, Length toUnit);
  static Uncertain<double> getFactor(Mass fromUnit, Mass toUnit);
  static Uncertain<double> getFactor(Viscosity fromUnit, Viscosity toUnit);

  static Uncertain<double> getOffset(Temperature fromUnit, Temperature toUnit);

  static std::string getName(Length unit);
  static std::string getName(Temperature unit);
  static std::string getName(Mass unit);
  static std::string getName(Viscosity unit);

  static Uncertain<double> kB();

private:
  Units() {}
};

// ================================================================

#endif  // #ifndef UNITS_H_

// ================================================================

// Local Variables:
// time-stamp-line-limit: 30
// mode: c++
// End:
