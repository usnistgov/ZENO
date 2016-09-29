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
// Date:    Mon Dec 29 11:44:20 2014 EDT
//
// Time-stamp: <2016-09-20 15:23:28 dcj>
//
// ================================================================

#include "SPRNG.h"

#include <sprng_cpp.h>
#include <cassert>

// ================================================================

/// Constructs a random number generator with the given stream number (out of
/// the given total number of streams) and the given seed.
///
SPRNG::SPRNG(int streamNum, int numStreams, int seed) {
  stream = SelectType(DEFAULT_RNG_TYPE);

  stream->init_sprng(streamNum, 
		     numStreams, 
		     seed, 
		     SPRNG_DEFAULT);

  assert(stream != NULL);
}

SPRNG::~SPRNG() {
  stream->free_sprng();
}

/// Return a random number in the range [0, 1).
///
double 
SPRNG::getRandIn01() {
  double num = stream->sprng();

  assert(num != -1);

  return num;
}

/// Return a random number in the range [min, max).
///
double 
SPRNG::getRandInRange(double min, double max) {
  double num = stream->sprng();

  assert(num != -1);

  return num * (max - min) + min;
}

// Local Variables:
// time-stamp-line-limit: 30
// mode: c++
// End:
