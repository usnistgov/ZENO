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

/// Defines a particle as an assembly of spheres with mutable center and orientation.
///

namespace zeno {

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

#include "Virials/Particle.h"

///  Constructs the class to define a particle as an assembly of spheres with mutable center and orientation.
///

using namespace zeno;

template <class T>
Particle<T>::
Particle(MixedModelProcessed<T> const & model, Sphere<T> const & boundingSphere) : model(model), boundingSphere(boundingSphere)
              {
}

template <class T>
Particle<T>::
  ~Particle() {
}

/// Obtain number of spheres in the particle.
///
template <class T>
int
Particle<T>::
numSpheres(){
    return model.getSpheres() -> size();
}

/// Obtain center of particle. Note: This is a constant, hence cannot be modified directly.
///
template <class T>
const Vector3<T>
Particle<T>::
getCenter() const {
    return center;
}

/// Sets the center of particle equal to the vector3 obtained as a function parameter.
///
template <class T>
void
Particle<T>::
setCenter(Vector3<T> v) {
    center = v;
}

/// Translates the particle by a step obtained as a function parameter.
///
template <class T>
void
Particle<T>::
translateBy(Vector3<T> step){
    center += step;
}

/// Rotates the particle about an axis and an by an angle, both obtained as a function parameters.
///
template <class T>
void
Particle<T>::
rotateBy(Vector3<T> axis, T angle){
    Matrix3x3<T> rotation;
    rotation.setAxisAngle(axis, angle);
    rotation.transform(orientation);
}

/// Sets the sphere at index of particle based on new position of particle.
///
template <class T>
const Vector3<T>
Particle<T>::
getSpherePosition( int index) const {
    Vector3<T> position = model.getSpheres() -> at(index).getCenter();
    orientation.transform(position);
    position += center;
    return position;
}

/// Obtain the assembly of spheres which constitute a particle.
///
template <class T>
MixedModelProcessed<T> const *
Particle<T>::
getModel(){
    return &model;
}

/// Returns a bounding sphere around the assembly of spheres.
///
template <class T>
Sphere<T> const *
Particle<T>::
getBoundingSphere(){
    return &boundingSphere;
}

template <class T>
const Vector3<T>
Particle<T>::
getBoundingSpherePosition() const {
    Vector3<T> position = boundingSphere.getCenter();
    orientation.transform(position);
    position += center;
    return position;
}

}

#endif
