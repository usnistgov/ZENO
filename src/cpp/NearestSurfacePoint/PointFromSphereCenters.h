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
// Date:    Mon Aug 17 17:01:47 2015 EDT
//
// Time-stamp: <2016-09-20 14:58:01 dcj>
//
// ================================================================

#ifndef POINT_FROM_SPHERE_CENTERS_H
#define POINT_FROM_SPHERE_CENTERS_H

#include "../Geometry/Sphere.h"
#include "../Geometry/Vector3.h"

// ================================================================

/// Finds the nearest point on the surface of an object represented by
/// sphere centers stored in spatial data structures, with spheres of each
/// radius being stored in a separate data structure.
///
template <class SphereCenterModel>
class PointFromSphereCenters
{
public:
  PointFromSphereCenters(SphereCenterModel const & sphereCenterModel);
  ~PointFromSphereCenters();

  void findNearestPoint(Vector3<double> const & queryPoint,
			double fracErrorBound,
			Vector3<double> * nearestPointNormal,
			double * distance) const;

private:
  SphereCenterModel const * sphereCenterModel;
};

template <class SphereCenterModel>
PointFromSphereCenters<SphereCenterModel>::
PointFromSphereCenters(SphereCenterModel const & sphereCenterModel) 
  : sphereCenterModel(&sphereCenterModel) {
  
}

template <class SphereCenterModel>
PointFromSphereCenters<SphereCenterModel>::
~PointFromSphereCenters() {

}

/// Compute the distance from the given query point to the nearest point on the 
/// surface of an object, allowing for the given relative error in distance.  
/// Also compute the surface normal at that point.
///
template <class SphereCenterModel>
void 
PointFromSphereCenters<SphereCenterModel>::
findNearestPoint(Vector3<double> const & queryPoint,
		 double fracErrorBound, 
		 Vector3<double> * nearestPointNormal,
		 double * distance) const {

  Sphere<double> const * nearestSphere = NULL;

  double minDistance = std::numeric_limits<double>::max();

  for (unsigned int radiusNum = 0; 
       radiusNum < sphereCenterModel->getNumRadii(); 
       radiusNum ++) {

    double centerDistSqr = -1;

    Sphere<double> const * foundSphere = NULL;

    sphereCenterModel->findNearestSphere(radiusNum,
					 queryPoint,
					 fracErrorBound,
					 &foundSphere, 
					 &centerDistSqr); 

    double foundDistance = 
      sqrt(centerDistSqr) - foundSphere->getRadius();

    if (foundDistance < minDistance) {
      nearestSphere = foundSphere;
      minDistance   = foundDistance;
    }
  }

  (*nearestPointNormal) = queryPoint - nearestSphere->getCenter();
  (*distance)           = minDistance;
}

// ================================================================

#endif  // #ifndef POINT_FROM_SPHERE_CENTERS_H

// ================================================================

// Local Variables:
// time-stamp-line-limit: 30
// mode: c++
// End:
