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

#ifndef RANDOM_UTILITIES_H
#define RANDOM_UTILITIES_H

#include "../Geometry/Vector3.h"

namespace zeno {

/// Generates random vectors inside and on a unit sphere.
///
template <class T, 
        class RandomNumberGenerator>
class RandomUtilities {
 public:
    RandomUtilities(RandomNumberGenerator * randomNumberGenerator);
  ~RandomUtilities();

  void setRandomOnSphere(Vector3<T> * v);
  void setRandomInSphere(Vector3<T> * v);

 private:
    RandomNumberGenerator * randomNumberGenerator;
};

}
#endif


