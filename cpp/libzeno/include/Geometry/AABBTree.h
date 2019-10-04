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
// Created: Fri Feb 10 11:41:46 2017 EDT
//
// ================================================================

#ifndef AABB_TREE_H
#define AABB_TREE_H

#include <vector>
#include <algorithm>

#include "Vector3.h"
#include "Cuboid.h"

// ================================================================

namespace zeno {
  
/// Spatial data structure for nearest-object queries based on an
/// Axis Aligned Bounding Box (AABB) tree.  The data structure can store any
/// objects which themselves can be enclosed in an AABB and which can compute
/// their distance to a point.  Specifically, the objects must provide the
/// following functions:
///
/// double getMaxCoord(int dim)
/// double getMinCoord(int dim)
/// double getDistanceSqrTo(Vector3<double> const & point)
/// bool contains(Vector3<double> const & point)
///
template <class ObjectT>
class AABBTree {  
 private:
  using objects_size_type  = typename std::vector<ObjectT>::size_type;
  
  struct AABBNode {
    Cuboid<double> aabb;

    objects_size_type objectsBegin;
    objects_size_type objectsEnd;
  };

  using aabbTree_size_type = typename std::vector<AABBNode>::size_type;

 public:
  AABBTree(objects_size_type maxNodeSize = 4);

  ~AABBTree();
  
  void preprocess(std::vector<ObjectT> * objects);

  void findNearestObject(Vector3<double> const & queryPoint,
			 ObjectT const * * nearestObject, 
			 double * minDistanceSqr) const;

  bool objectsContain(Vector3<double> const & queryPoint) const;

  bool isEmpty() const;
  
  void printDataStructureStats() const;
  void printSearchStats() const;

 private:
  const objects_size_type maxNodeSize;
  
  std::vector<AABBNode> aabbTree;

  std::vector<ObjectT> * sortedObjects;

  void computeAABBBounds(aabbTree_size_type nodeIndex);

  void splitAABBNode(aabbTree_size_type nodeIndex,
		     int splitDim);
};

// ================================================================

/// Constructs an empty spatial data structure.
///
template <class ObjectT>
AABBTree<ObjectT>::AABBTree(objects_size_type maxNodeSize)
  : maxNodeSize(maxNodeSize),
    aabbTree(),
    sortedObjects(NULL) {

}

template <class ObjectT>
AABBTree<ObjectT>::~AABBTree() {

}

/// Inserts the given objects into the spatial data structure.
/// The objects vector may be modified and must not be changed outside the
/// class.
///
template <class ObjectT>
void
AABBTree<ObjectT>::preprocess(std::vector<ObjectT> * objects) {
  
  aabbTree.clear();
  
  sortedObjects = objects;

  if (objects->empty()) {
    return;
  }
  
  aabbTree.resize(1);

  aabbTree.at(0).objectsBegin = 0;
  aabbTree.at(0).objectsEnd   = sortedObjects->size();

  aabbTree_size_type nodesBegin = 0;
  aabbTree_size_type nodesEnd   = 1;

  int splitDim = 0;

  bool splitOccurred = true;
  
  while (splitOccurred) {
    splitOccurred = false;
    
    for (aabbTree_size_type nodeIndex = nodesBegin;
	 nodeIndex < nodesEnd;
	 ++nodeIndex) {

      computeAABBBounds(nodeIndex);

      objects_size_type nodeSize =
	aabbTree.at(nodeIndex).objectsEnd -
	aabbTree.at(nodeIndex).objectsBegin;
      
      if (nodeSize > maxNodeSize) {
	if (!splitOccurred) {
	  aabbTree_size_type newAABBTreeSize = aabbTree.size() * 2 + 1; 

	  aabbTree.resize(newAABBTreeSize);
	}

	splitAABBNode(nodeIndex, splitDim);
	
	splitOccurred = true;
      }
    }

    nodesBegin = nodesEnd;
    nodesEnd   = nodesEnd * 2 + 1;

    ++ splitDim;

    if (splitDim == 3) {
      splitDim = 0;
    }
  }
}

template <class ObjectT>
void
AABBTree<ObjectT>::computeAABBBounds(aabbTree_size_type nodeIndex) {
  for (int dim = 0; dim < 3; ++dim) {
    aabbTree.at(nodeIndex).aabb.
      setMaxCoord(dim, std::numeric_limits<double>::lowest());
  }

  for (int dim = 0; dim < 3; ++dim) {
    aabbTree.at(nodeIndex).aabb.
      setMinCoord(dim, std::numeric_limits<double>::max());
  }

  for (objects_size_type objectIndex =
	 aabbTree.at(nodeIndex).objectsBegin;
       objectIndex < aabbTree.at(nodeIndex).objectsEnd;
       ++ objectIndex) {

    for (int dim = 0; dim < 3; ++dim) {
      aabbTree.at(nodeIndex).aabb.
	setMaxCoord(dim,
		    std::max(aabbTree.at(nodeIndex).aabb.getMaxCoord(dim),
			     sortedObjects->at(objectIndex).getMaxCoord(dim)));
    }

    for (int dim = 0; dim < 3; ++dim) {
      aabbTree.at(nodeIndex).aabb.
	setMinCoord(dim,
		    std::min(aabbTree.at(nodeIndex).aabb.getMinCoord(dim),
			     sortedObjects->at(objectIndex).getMinCoord(dim)));
    }
  }
}

template <class ObjectT>
void
AABBTree<ObjectT>::splitAABBNode(aabbTree_size_type nodeIndex,
				 int splitDim) {

  class CoordComp {
   public:
    CoordComp(int dim) : dim(dim) {}
    
    bool operator() (ObjectT i, ObjectT j) {
      double iMidCoordX2 = i.getMaxCoord(dim) + i.getMinCoord(dim);
      double jMidCoordX2 = j.getMaxCoord(dim) + j.getMinCoord(dim);

      return (iMidCoordX2 < jMidCoordX2);
    }

   private:
    int dim; 
  };
  
  CoordComp coordComp(splitDim);
	
  std::sort(sortedObjects->begin() + aabbTree.at(nodeIndex).objectsBegin,
	    sortedObjects->begin() + aabbTree.at(nodeIndex).objectsEnd,
	    coordComp);

  objects_size_type nodeSize =
    aabbTree.at(nodeIndex).objectsEnd -
    aabbTree.at(nodeIndex).objectsBegin;
  
  objects_size_type objectsMiddle =
    aabbTree.at(nodeIndex).objectsBegin +
    nodeSize / 2;

  aabbTree_size_type lowChildIndex  = nodeIndex * 2 + 1;
  aabbTree_size_type highChildIndex = nodeIndex * 2 + 2;

  aabbTree.at(lowChildIndex).objectsBegin =
    aabbTree.at(nodeIndex).objectsBegin;

  aabbTree.at(lowChildIndex).objectsEnd =
    objectsMiddle;

  aabbTree.at(highChildIndex).objectsBegin =
    objectsMiddle;

  aabbTree.at(highChildIndex).objectsEnd =
    aabbTree.at(nodeIndex).objectsEnd;

  aabbTree.at(nodeIndex).objectsBegin = 0;
  aabbTree.at(nodeIndex).objectsEnd   = 0;
}

/// Searches the objects for the nearest to the given query point.
/// Computes both the nearest object and its distance squared.
/// If no objects are found that are nearer than the initial minDistanceSqr,
/// nearestObject and minDistanceSqr are not modified.
///
template <class ObjectT>
void
AABBTree<ObjectT>::findNearestObject(Vector3<double> const & queryPoint,
				     ObjectT const * * nearestObject, 
				     double * minDistanceSqr) const {

  if (aabbTree.empty()) {
    return;
  }
  
  struct NodeToCheck {
    NodeToCheck(aabbTree_size_type nodeIndex,
		double distanceSqrToNode)
      : nodeIndex(nodeIndex),
	distanceSqrToNode(distanceSqrToNode) {}
    
    aabbTree_size_type nodeIndex;
    
    double distanceSqrToNode;
  };

  std::vector<NodeToCheck> nodesToCheck;

  nodesToCheck.emplace_back(0,
			    aabbTree.at(0).aabb.getDistanceSqrTo(queryPoint));
  
  while(!nodesToCheck.empty()) {
    aabbTree_size_type nodeIndex = nodesToCheck.back().nodeIndex;
    double distanceSqrToNode     = nodesToCheck.back().distanceSqrToNode;

    nodesToCheck.pop_back();

    if (*minDistanceSqr < distanceSqrToNode) {
      continue;
    }

    objects_size_type objectsBegin = aabbTree.at(nodeIndex).objectsBegin;
    objects_size_type objectsEnd   = aabbTree.at(nodeIndex).objectsEnd;

    if (objectsBegin == objectsEnd) {
      // interior node

      aabbTree_size_type lowChildIndex  = nodeIndex * 2 + 1;
      aabbTree_size_type highChildIndex = nodeIndex * 2 + 2;

      double distanceSqrToLowChild =
	aabbTree.at(lowChildIndex).aabb.getDistanceSqrTo(queryPoint);

      double distanceSqrToHighChild =
	aabbTree.at(highChildIndex).aabb.getDistanceSqrTo(queryPoint);

      if (distanceSqrToLowChild < distanceSqrToHighChild) {
	// check low child first
	nodesToCheck.emplace_back(highChildIndex, distanceSqrToHighChild);
	nodesToCheck.emplace_back(lowChildIndex, distanceSqrToLowChild);
      }
      else {
	// check high child first
	nodesToCheck.emplace_back(lowChildIndex, distanceSqrToLowChild);
	nodesToCheck.emplace_back(highChildIndex, distanceSqrToHighChild);
      }
    }
    else {
      // leaf node

      for (objects_size_type objectIndex = objectsBegin;
	   objectIndex < objectsEnd;
	   ++objectIndex) {

	double distanceSqrToObject =
	  sortedObjects->at(objectIndex).getDistanceSqrTo(queryPoint);

	if (distanceSqrToObject < *minDistanceSqr) {
	  *nearestObject  = &(sortedObjects->at(objectIndex));
	  *minDistanceSqr = distanceSqrToObject;
	}
      }
    }
  }
}

/// Returns whether any object in the tree contains the given query point.
///
template <class ObjectT>
bool
AABBTree<ObjectT>::objectsContain(Vector3<double> const & queryPoint)
  const {

  if (aabbTree.empty()) {
    return false;
  }
  
  std::vector<aabbTree_size_type> nodesToCheck;

  nodesToCheck.emplace_back(0);

  while(!nodesToCheck.empty()) {
    aabbTree_size_type nodeIndex = nodesToCheck.back();

    nodesToCheck.pop_back();

    objects_size_type objectsBegin = aabbTree.at(nodeIndex).objectsBegin;
    objects_size_type objectsEnd   = aabbTree.at(nodeIndex).objectsEnd;

    if (objectsBegin == objectsEnd) {
      // interior node

      aabbTree_size_type lowChildIndex  = nodeIndex * 2 + 1;
      aabbTree_size_type highChildIndex = nodeIndex * 2 + 2;

      if (aabbTree.at(highChildIndex).aabb.contains(queryPoint)) {
	nodesToCheck.emplace_back(highChildIndex);
      }

      if (aabbTree.at(lowChildIndex).aabb.contains(queryPoint)) {
	nodesToCheck.emplace_back(lowChildIndex);
      }
    }
    else {
      // leaf node

      for (objects_size_type objectIndex = objectsBegin;
	   objectIndex < objectsEnd;
	   ++objectIndex) {

	if (sortedObjects->at(objectIndex).contains(queryPoint)) {
	  return true;
	}
      }
    }
  }

  return false;
}

template <class ObjectT>
bool
AABBTree<ObjectT>::isEmpty() const {
  return aabbTree.empty();
}

template <class ObjectT>
void
AABBTree<ObjectT>::printDataStructureStats() const {

}

template <class ObjectT>
void
AABBTree<ObjectT>::printSearchStats() const {

}

}

#endif
