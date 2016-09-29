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
// Date:    Wed May 21 16:21:57 2014 EDT
//
// Time-stamp: <2016-08-30 11:45:53 dcj>
//
// ================================================================

#ifndef TIMER_H_
#define TIMER_H_

// ================================================================

#include <chrono>

// ================================================================

/// A stopwatch for benchmarking.
///
class Timer
{
public:
  Timer();
  ~Timer();

  void start();
  void stop();
  void reset();

  double getTime() const;

private:
  bool running;

  std::chrono::high_resolution_clock::time_point startTime;

  double elapsedTime;
};
// ================================================================

#endif  // #ifndef TIMER_H_

// ================================================================

// Local Variables:
// time-stamp-line-limit: 30
// mode: c++
// End:
