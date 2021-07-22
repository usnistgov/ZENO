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

#ifndef CLUSTER_SUM_H
#define CLUSTER_SUM_H

#include "OverlapTester.h"

namespace zeno {

/// Computes cluster sum.
///

template <class T>
class ClusterSum {
 public:
    ClusterSum(std::vector<Particle<T> *> * particles, OverlapTester<T> const * overlapTester);
    virtual ~ClusterSum();
    virtual double value() = 0;

 protected:
    std::vector<Particle<T> *> * particles;
    OverlapTester<T> const * overlapTester;
};

///Sub class of ClusterSum to compute cluster sum for chains.
///

template <class T>
class ClusterSumChain : public ClusterSum<T> {
public:
    ClusterSumChain(std::vector<Particle<T> *> * particles, double diameter, double ringFac, double chainFac);
    ~ClusterSumChain();
    double value();

private:
      double diameter;
      double ringFac;
      double chainFac;
};

///Sub class of ClusterSum to compute cluster sum using Wheatley Recursion.
///
template <class T>
class ClusterSumWheatleyRecursion : public ClusterSum<T>{
public:
    ClusterSumWheatleyRecursion(std::vector<Particle<T> *> * particles, OverlapTester<T> const * overlapTester);
    ~ClusterSumWheatleyRecursion();
    double value();

private:
    double preFac;
};

}
#endif

