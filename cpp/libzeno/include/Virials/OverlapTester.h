// ================================================================
//
/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
// 
// ================================================================

// ================================================================
// 
// Authors:
// Created:
//
// ================================================================

#ifndef OVERLAP_TESTER_H
#define OVERLAP_TESTER_H

#include "Particle.h"

namespace zeno {

/// Checks if two particles are overlapped.
///
template <class T>
class OverlapTester {
 public:
    OverlapTester();
  ~OverlapTester();

  bool isOverlapped(Particle<T> * a, Particle<T> * b) const;
};

}
#endif

