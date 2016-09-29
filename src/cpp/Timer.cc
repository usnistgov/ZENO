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
// Date:    Wed May 21 16:22:30 2014 EDT
// 
// Time-stamp: <2016-08-30 11:50:34 dcj>
// 
// ================================================================

#include "Timer.h"

#include <ctime>
#include <ratio>

// ================================================================

using namespace std::chrono;

// ================================================================

/// Construct a Timer with 0s on the clock.
///
Timer::Timer() 
  : running(false),
    startTime(),
    elapsedTime(0) {

}

Timer::~Timer() {

}

/// Start the Timer.  Any previous time on the Timer is not cleared.
///
void 
Timer::start() {
  if (!running) {
    running = true;

    startTime = high_resolution_clock::now();
  }
}

/// Stop the Timer.  The time on the Timer is not cleared.
///
void 
Timer::stop() {
  if (running) {
    running = false;

    high_resolution_clock::time_point endTime = high_resolution_clock::now();

    elapsedTime += 
      (duration_cast<duration<double>>(endTime - startTime)).count();
  }
}

/// Set the time on the Timer to 0s.  If the Timer is running, it is not
/// stopped.
///
void 
Timer::reset() {
  elapsedTime = 0;

  startTime = high_resolution_clock::now();
}


/// Returns the current time on the Timer.  Can be called even if the Timer is
/// running.
///
double
Timer::getTime() const {
  if (running) {
    high_resolution_clock::time_point endTime = high_resolution_clock::now();

    double currentTime = 
      (duration_cast<duration<double>>(endTime - startTime)).count();

    return currentTime + elapsedTime;
  }
  else {
    return elapsedTime;
  }
}

// ================================================================

// Local Variables:
// time-stamp-line-limit: 30
// End:
