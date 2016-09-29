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
// Date:    Wed Nov 05 15:29:00 2014 EDT
//
// Time-stamp: <2016-09-20 16:19:40 dcj>
//
// ================================================================

#include <vector>
#include <cassert>

#include <nanoflann.hpp>

#include "../Geometry/Sphere.h"
#include "../Geometry/Vector3.h"

// ================================================================

#ifndef NANOFLANNDATASETADAPTOR_H
#define NANOFLANNDATASETADAPTOR_H

/// Adaptor class used by the nanoFLANN library to access sphere center data.
///
template <typename DistanceType, typename ComponentType>
class NanoFLANNDatasetAdaptor {
 public:
  NanoFLANNDatasetAdaptor(int nPts,
			  int * originalSphereIndexes,
			  std::vector<Sphere<double> > const * originalSpheres);

  ~NanoFLANNDatasetAdaptor();

  /// Must return the number of data points
  size_t kdtree_get_point_count() const;

  /// Must return the Euclidean (L2) distance between the vector "p1[0:size-1]"
  /// and the data point with index "idx_p2" stored in the class:
  DistanceType kdtree_distance(const ComponentType *p1, const size_t idx_p2, 
			       size_t size) const;

  /// Must return the dim'th component of the idx'th point in the class:
  ComponentType kdtree_get_pt(const size_t idx, int dim) const;

  /// Optional bounding-box computation: return false to default to a standard 
  /// bbox computation loop.
  /// Return true if the BBOX was already computed by the class and returned in
  /// "bb" so it can be avoided to redo it again.
  /// Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 
  /// for point clouds)
  template <class BBOX>
  bool kdtree_get_bbox(BBOX &bb) const;

 private:
  int nPts;

  ComponentType * points;
};

template <typename DistanceType, typename ComponentType>
NanoFLANNDatasetAdaptor<DistanceType, ComponentType>::
NanoFLANNDatasetAdaptor(int nPts,
			int * originalSphereIndexes, 
			std::vector<Sphere<double> > const * originalSpheres)
  : nPts(nPts),
    points(NULL) {

  points = new ComponentType[nPts * 3];

  assert(points != NULL);

  for (int i = 0; i < nPts; i++) {
    const Vector3<double> center = 
      originalSpheres->at(originalSphereIndexes[i]).getCenter();

    for (int dim = 0; dim < 3; dim++) {
      points[i*3 + dim] = center.get(dim);
    }
  }
}

template <typename DistanceType, typename ComponentType>
NanoFLANNDatasetAdaptor<DistanceType, ComponentType>::
  ~NanoFLANNDatasetAdaptor () {

  delete [] points;
}

template <typename DistanceType, typename ComponentType>
size_t 
NanoFLANNDatasetAdaptor<DistanceType, ComponentType>::
  kdtree_get_point_count() const {

  return nPts;
}

template <typename DistanceType, typename ComponentType>
DistanceType 
NanoFLANNDatasetAdaptor<DistanceType, ComponentType>::
  kdtree_distance(const ComponentType *p1, 
		  const size_t idx_p2, 
		  size_t size) const {

  assert(size == 3);

  DistanceType distSqr = 0;

  for (int dim = 0; dim < 3; dim++) {
    distSqr += pow(p1[dim] - points[idx_p2*3 + dim], 2);
  }

  return distSqr;
}

template <typename DistanceType, typename ComponentType>
ComponentType 
NanoFLANNDatasetAdaptor<DistanceType, ComponentType>::
  kdtree_get_pt(const size_t idx, 
		int dim) const {

  return points[idx*3 + dim];
}

template <typename DistanceType, typename ComponentType>
template <class BBOX>
bool 
NanoFLANNDatasetAdaptor<DistanceType, ComponentType>::
  kdtree_get_bbox(BBOX &bb) const {

  return false;
}

#endif

// ================================================================

// Local Variables:
// time-stamp-line-limit: 30
// End:
