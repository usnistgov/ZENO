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
// Created: Wed May 21 16:22:30 2014 EDT
// 
// ================================================================

#include "Timer.h"

#include <ctime>
#include <ratio>

// ================================================================

using namespace std::chrono;

using namespace zeno;

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


/// Returns the current time on the Timer (in seconds).  Can be called even if
/// the Timer is running.
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

