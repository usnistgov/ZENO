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
// Created: 2019-07-03
// 
// ================================================================

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "Virials/ParametersVirial.h"

// ================================================================

using namespace zeno;

ParametersVirial::ParametersVirial() 
  : numThreads(),
    seed(),
    steps(),
    stepsWasSet(false),
    order(),
    orderWasSet(false) {

}

ParametersVirial::~ParametersVirial() {

}

void
ParametersVirial::setNumThreads(int numThreads) {
  this->numThreads = numThreads;
}

int 
ParametersVirial::getNumThreads() const {
  return numThreads;
}

void
ParametersVirial::setSeed(int seed) {
  this->seed = seed;
}

int 
ParametersVirial::getSeed() const {
  return seed;
}

void
ParametersVirial::setOrder(int order) {
  this->order = order;

  orderWasSet = true;
}

int
ParametersVirial::getOrder() const {
  return order;
}

bool
ParametersVirial::getOrderWasSet() const {
  return orderWasSet;
}

void
ParametersVirial::setSteps(long long steps) {
  this->steps = steps;

  stepsWasSet = true;
}

long long 
ParametersVirial::getSteps() const {
  return steps;
}

bool 
ParametersVirial::getStepsWasSet() const {
  return stepsWasSet;
}

void 
ParametersVirial::mpiBroadcast(int root) {
#ifdef USE_MPI
  int mpiSize = 1, mpiRank = 0;

  MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

  if (mpiSize > 1) {
    if (mpiRank == root) {
      serializeMpiBroadcast(root);
    }
    else {
      mpiBroadcastDeserialize(root);
    }
  }
#endif
}

/// Broadcasts the parameters that cannot be set in the bod file over MPI.
///
void 
ParametersVirial::serializeMpiBroadcast(int root) const {
#ifdef USE_MPI
  const int numLongLongsToSend = 5;

  long long longLongsArray[numLongLongsToSend];

  longLongsArray[0] = (long long)getNumThreads();
  longLongsArray[1] = (long long)getSeed();
  longLongsArray[2] = (long long)getStepsWasSet();
  longLongsArray[3] = getSteps();
  longLongsArray[4] = getOrder();

  MPI_Bcast(longLongsArray, numLongLongsToSend, MPI_LONG_LONG,
	    root, MPI_COMM_WORLD);
#endif
}

/// Receives the parameters that cannot be set in the bod file over MPI.
///
void 
ParametersVirial::mpiBroadcastDeserialize(int root) {
#ifdef USE_MPI
  const int numLongLongsToReceive = 5;

  long long longLongsArray[numLongLongsToReceive];

  MPI_Bcast(longLongsArray, numLongLongsToReceive, MPI_LONG_LONG,
	    root, MPI_COMM_WORLD);

  setNumThreads((int)longLongsArray[0]);

  setSeed((int)longLongsArray[1]);

  if ((bool)longLongsArray[3]) {
    setSteps(longLongsArray[2]);
  }

  setOrder((bool)longLongsArray[4]);
#endif
}
