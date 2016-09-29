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
// Date:    Fri Aug 15 14:35:03 2014 EDT
//
// Time-stamp: <2016-09-20 16:48:15 dcj>
//
// ================================================================

#ifndef NANOFLANNSORT_H
#define NANOFLANNSORT_H

#include <vector>

#include <nanoflann.hpp>

#include "NanoFLANNDatasetAdaptor.h"

#include "../Geometry/Sphere.h"
#include "../Geometry/Vector3.h"

// ================================================================

/// Spatial data structure based on the nanoFLANN library for storing spheres.
/// Spheres are sorted by radius, and the centers of the spheres of each radius
/// are stored in a separate nanoFLANN data structure.
///
class NanoFLANNSort {
 public:
  NanoFLANNSort();

  ~NanoFLANNSort();

  void preprocess(std::vector<Sphere<double> > const & spheres,
		  double fracErrorBound);

  unsigned int getNumRadii() const;

  void findNearestSphere(int radiusNum,
			 Vector3<double> const & queryPoint,
			 double fracErrorBound,
			 Sphere<double> const * * nearestSphere, 
			 double * centerDistSqr) const;

  void printDataStructureStats() const;
  void printSearchStats() const;

 private:
  typedef NanoFLANNDatasetAdaptor<double, double> DatasetAdaptorType;

  typedef nanoflann::L2_Simple_Adaptor<double, DatasetAdaptorType> MetricType;

  typedef nanoflann::KDTreeSingleIndexAdaptor<MetricType, 
                                              DatasetAdaptorType, 3, int> 
    KDTreeType;

  struct NanoFLANNInstance {
    double sphereRadiusSqr;

    int nPts;
    
    DatasetAdaptorType * dataset;
    KDTreeType * kdTree;

    int * originalSphereIndexes;
  };

  const std::vector<Sphere<double> > * originalSpheres;

  std::vector<NanoFLANNInstance> nanoFLANNInstances;
};

#endif

// ================================================================

// Local Variables:
// time-stamp-line-limit: 30
// End:
