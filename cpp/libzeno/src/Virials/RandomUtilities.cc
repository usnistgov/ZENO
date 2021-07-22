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

#include <cmath>

#include "RandomUtilities.h"

/// Constructs the class to generate random vectors inside and on a unit sphere.
///
template <class T,
        class RandomNumberGenerator>
RandomUtilities<T, RandomNumberGenerator>::
RandomUtilities(RandomNumberGenerator * randomNumberGenerator):randomNumberGenerator(randomNumberGenerator){
}

template <class T,
        class RandomNumberGenerator>
RandomUtilities<T, RandomNumberGenerator>::
  ~RandomUtilities() {
}

/// Generates random vectors on a unit sphere.
///
template <class T,
        class RandomNumberGenerator>
void
RandomUtilities<T, RandomNumberGenerator>::
setRandomOnSphere(Vector3<T> * v) {
    double z1, z2, zsq;
    do  {
        z1 = 2.0 * randomNumberGenerator->getRandIn01() - 1.0;
        z2 = 2.0 * randomNumberGenerator->getRandIn01() - 1.0;
        zsq = z1 * z1 + z2 * z2;
    } while (zsq > 1.0);

    double ranh = 2.0 * sqrt(1.0 - zsq);
    v->setXYZ(z1 * ranh, z2 * ranh, 1.0 - 2.0 * zsq);
}

/// Generates random vectors inside a unit sphere.
///
template <class T,
        class RandomNumberGenerator>
void
RandomUtilities<T, RandomNumberGenerator>::
setRandomInSphere(Vector3<T> * v) {
    double r = cbrt(randomNumberGenerator->getRandIn01());
    double u, w, s;
    do {
        u = 1.0 - 2.0*randomNumberGenerator->getRandIn01();
        w = 1.0 - 2.0*randomNumberGenerator->getRandIn01();
        s = u*u + w*w;
    } while(s > 1);

    double ra = 2 * r * sqrt(1 - s);
    v->setXYZ(ra * u, ra * w, r * (2 * s - 1));

}

