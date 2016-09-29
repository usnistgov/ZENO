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
// Date:    Fri Aug 15 14:35:07 2014 EDT
//
// Time-stamp: <2016-09-20 17:01:57 dcj>
//
// ================================================================

#include <limits>
#include <cassert>
#include <iostream>

#include "NanoFLANNSort.h"

// ================================================================

/// Constructs an empty spatial data structure.
///
NanoFLANNSort::NanoFLANNSort()
  : originalSpheres(NULL),
    nanoFLANNInstances() {

}

NanoFLANNSort::~NanoFLANNSort() {
  for (unsigned int i = 0; i < nanoFLANNInstances.size(); i ++) {
    delete nanoFLANNInstances.at(i).dataset;

    delete nanoFLANNInstances.at(i).kdTree;

    delete [] nanoFLANNInstances.at(i).originalSphereIndexes;
  }
}

/// Inserts the given spheres into the spatial data structure, allowing for
/// the given relative error in distance.
///
void 
NanoFLANNSort::preprocess(std::vector<Sphere<double> > const & spheres,
			  double fracErrorBound) {

  originalSpheres = &spheres;

  //count number of spheres of each radius
  for (unsigned int sphereNum = 0; 
       sphereNum < spheres.size(); 
       sphereNum ++) {

    double sphereRadiusSqr = spheres.at(sphereNum).getRadiusSqr();

    bool radiusFound = false;

    for (unsigned int radiusNum = 0; 
	 radiusNum < nanoFLANNInstances.size(); 
	 radiusNum ++) {

      if (nanoFLANNInstances.at(radiusNum).sphereRadiusSqr == sphereRadiusSqr) {
	radiusFound = true;

	nanoFLANNInstances.at(radiusNum).nPts ++;

	break;
      }
    }

    if (!radiusFound) {
      NanoFLANNInstance nanoFLANNInstance;

      nanoFLANNInstance.sphereRadiusSqr = sphereRadiusSqr;
    
      nanoFLANNInstance.nPts = 1;

      nanoFLANNInstance.dataset = NULL;
      nanoFLANNInstance.kdTree  = NULL;

      nanoFLANNInstance.originalSphereIndexes = NULL;

      nanoFLANNInstances.push_back(nanoFLANNInstance);
    }
  }

  //allocate memory for points
  for (unsigned int radiusNum = 0; 
       radiusNum < nanoFLANNInstances.size(); 
       radiusNum ++) {

    int nPts = nanoFLANNInstances.at(radiusNum).nPts;

    nanoFLANNInstances.at(radiusNum).nPts = 0; //reset point counter

    nanoFLANNInstances.at(radiusNum).originalSphereIndexes = new int[nPts];
  }

  //add spheres to instances
  for (unsigned int sphereNum = 0; 
       sphereNum < spheres.size(); 
       sphereNum ++) {

    double sphereRadiusSqr = spheres.at(sphereNum).getRadiusSqr();

    bool radiusFound = false;

    for (unsigned int radiusNum = 0; 
	 radiusNum < nanoFLANNInstances.size(); 
	 radiusNum ++) {

      if (nanoFLANNInstances.at(radiusNum).sphereRadiusSqr == sphereRadiusSqr) {
	radiusFound = true;

	int nPts = nanoFLANNInstances.at(radiusNum).nPts;

	nanoFLANNInstances.at(radiusNum).originalSphereIndexes[nPts] = 
	  sphereNum;

	nanoFLANNInstances.at(radiusNum).nPts ++;

	break;
      }
    }

    assert(radiusFound);
  }

  //build kdTrees
  for (unsigned int radiusNum = 0; 
       radiusNum < nanoFLANNInstances.size(); 
       radiusNum ++) {

    int nPts = nanoFLANNInstances.at(radiusNum).nPts;

    int * originalSphereIndexes = 
      nanoFLANNInstances.at(radiusNum).originalSphereIndexes;

    DatasetAdaptorType * dataset = 
      new DatasetAdaptorType(nPts, originalSphereIndexes, originalSpheres);

    nanoFLANNInstances.at(radiusNum).dataset = dataset;

    nanoFLANNInstances.at(radiusNum).kdTree = new KDTreeType(3, *dataset);

    nanoFLANNInstances.at(radiusNum).kdTree->buildIndex();
  }
}

void 
NanoFLANNSort::printDataStructureStats() 
  const {

  std::cout << "Number of NanoFLANN instances: " << nanoFLANNInstances.size()
	    << std::endl
	    << "Used memory (MB)" << std::endl;

  for (unsigned int radiusNum = 0; 
       radiusNum < nanoFLANNInstances.size(); 
       radiusNum ++) {

    std::cout << radiusNum << ": " 
	      << nanoFLANNInstances.at(radiusNum).kdTree->usedMemory() / 1000000.
	      << std::endl;
  }

  std::cout << std::endl;
}

void 
NanoFLANNSort::printSearchStats() 
  const {

}

/// Return the number of unique sphere radii in the spatial data structure.
///
unsigned int 
NanoFLANNSort::
getNumRadii() 
  const {

  return nanoFLANNInstances.size();
}

/// Searches the spheres in the given radius bin for the closest to the given
/// query point, allowing for the given relative error in distance.  Computes
/// both the closest sphere and the square distance to its center.
///
void 
NanoFLANNSort::
findNearestSphere(int radiusNum,
		  Vector3<double> const & queryPoint,
		  double fracErrorBound,
		  Sphere<double> const * * nearestSphere, 
		  double * centerDistSqr) 
  const {

  double query[3];

  for (int i = 0; i < 3; i++) {
    query[i] = queryPoint.get(i);
  }

  int resultIndex = -1;

  // nanoFLANNInstances.at(radiusNum).kdTree->knnSearch(query, 1, 
  // 						&resultIndex, &distSqr);

  nanoflann::SearchParams params;
  params.eps = fracErrorBound;

  nanoflann::KNNResultSet<double, int> resultSet(1);
  resultSet.init(&resultIndex, centerDistSqr);

  nanoFLANNInstances.at(radiusNum).kdTree->findNeighbors(resultSet, query, 
						     params);

  int sphereIndex = 
    nanoFLANNInstances.at(radiusNum).originalSphereIndexes[resultIndex];

  *nearestSphere = &(originalSpheres->at(sphereIndex));
}

// ================================================================

// Local Variables:
// time-stamp-line-limit: 30
// End:
