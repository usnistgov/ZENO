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
// Time-stamp: <2016-09-20 15:02:25 dcj>
//
// ================================================================

#ifndef SPRNG_H_
#define SPRNG_H_

class Sprng;

// ================================================================

/// Generates random numbers using the SPRNG library.
///
class SPRNG
{
public:
  SPRNG(int streamNum, int numStreams, int seed);
  ~SPRNG();

  double getRandIn01();
  double getRandInRange(double min, double max);

private:
  Sprng *stream;
};

// ================================================================

#endif  // #ifndef SPRNG_H_

// ================================================================

// Local Variables:
// time-stamp-line-limit: 30
// mode: c++
// End:
