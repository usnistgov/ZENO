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

#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>

#include "../Geometry/MixedModelProcessed.h"
#include "../Geometry/Sphere.h"
#include "../Geometry/Vector3.h"
#include "../Geometry/Matrix3x3.h"

namespace zeno {

/// Defines a particle as an assembly of spheres with mutable center and orientation.
///
template <class T>
class Particle {
 public:
    Particle(MixedModelProcessed<T> const & model,
            Sphere<T> const & boundingSphere);
  ~Particle();

  int numSpheres();
  const Vector3<T> getCenter() const;
  void setCenter(Vector3<T> v);
  void translateBy(Vector3<T> step);
  void rotateBy(Vector3<T> axis, T angle);
  const Vector3<T> getSpherePosition(int index) const;
  const Vector3<T> getBoundingSpherePosition() const;
  MixedModelProcessed<T> const * getModel();
  Sphere<T> const * getBoundingSphere();

 private:
    MixedModelProcessed<T> const & model;
    Matrix3x3<T> orientation{1,0,0,0,1,0,0,0,1};
    Vector3<T> center{0,0,0};
    Sphere<T> const & boundingSphere;
};

}
#endif

