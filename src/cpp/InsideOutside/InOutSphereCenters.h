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
// Time-stamp: <2016-09-20 14:36:04 dcj>
//
// ================================================================

#ifndef IN_OUT_SPHERE_CENTERS_H
#define IN_OUT_SPHERE_CENTERS_H

#include "../Geometry/Vector3.h"

// ================================================================

/// Tests whether a point is inside or outside an object represented by
/// sphere centers stored in spatial data structures, with spheres of each
/// radius being stored in a separate data structure.
///
template <class SphereCenterModel>
class InOutSphereCenters
{
public:
  InOutSphereCenters(SphereCenterModel const & sphereCenterModel);
  ~InOutSphereCenters();

  bool isInside(Vector3<double> const & queryPoint,
		double fracErrorBound) const;

private:
  SphereCenterModel const * sphereCenterModel;
};

template <class SphereCenterModel>
InOutSphereCenters<SphereCenterModel>::
InOutSphereCenters(SphereCenterModel const & sphereCenterModel) 
  : sphereCenterModel(&sphereCenterModel) {
  
}

template <class SphereCenterModel>
InOutSphereCenters<SphereCenterModel>::
~InOutSphereCenters() {

}

/// Return whether the given query point is inside the object, allowing for
/// the given relative error in distance.
///
template <class SphereCenterModel>
bool 
InOutSphereCenters<SphereCenterModel>::
isInside(Vector3<double> const & queryPoint,
	 double fracErrorBound) const {

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

    if (centerDistSqr < foundSphere->getRadiusSqr()) {
      return true;
    }
  }

  return false;
}

// ================================================================

#endif  // #ifndef IN_OUT_SPHERE_CENTERS_H

// ================================================================

// Local Variables:
// time-stamp-line-limit: 30
// mode: c++
// End:
