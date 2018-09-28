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
// Created: Thu Feb 12 17:24:57 2015 EDT
//
// ================================================================

#ifndef QUASI_RANDOM_SPHERE_POINTS_REGULAR_CUBE_H
#define QUASI_RANDOM_SPHERE_POINTS_REGULAR_CUBE_H

#include <vector>
#include <cmath>

#include "../Geometry/Sphere.h"
#include "../Geometry/Vector3.h"

template <class T>
class QuasiRandomSpherePointsRegularCube {
 public:
  static void generate(Sphere<T> const & sphere,
		       int minNumPoints,
		       std::vector<Vector3<T> > * points);
};

template <class T>
void
QuasiRandomSpherePointsRegularCube<T>::
  generate(Sphere<T> const & sphere,
	   int minNumPoints,
	   std::vector<Vector3<T> > * points) {

  points->clear();

  //generate points on unit cube centered at origin

  double pointsPerFace = minNumPoints / 6.;
  double pointsPerEdge = ceil(sqrt(pointsPerFace));

  double pointSpacing = 1 / pointsPerEdge;

  for (double u = (pointSpacing / 2) - 0.5; u < 0.5; u += pointSpacing) {
    for (double v = (pointSpacing / 2) - 0.5; v < 0.5; v += pointSpacing) {

      points->push_back(Vector3<T>(u, v, -0.5));
      points->push_back(Vector3<T>(u, v,  0.5));
      points->push_back(Vector3<T>(u, -0.5, v));
      points->push_back(Vector3<T>(u,  0.5, v));
      points->push_back(Vector3<T>(-0.5, u, v));
      points->push_back(Vector3<T>( 0.5, u, v));
    }
  }

  for (typename std::vector<Vector3<T> >::iterator pointIt = points->begin(); 
       pointIt != points->end();
       ++pointIt) {

    //project point onto unit sphere centered at origin 

    pointIt->normalize();

    //move point onto requested sphere

    (*pointIt) *= sphere.getRadius();
    (*pointIt) += sphere.getCenter();
  }
}

#endif

