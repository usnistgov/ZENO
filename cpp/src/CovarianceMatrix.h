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
// Created: Wed Feb 17 09:02:26 2016 EDT
//
// ================================================================

#ifndef COVARIANCE_MATRIX_H_
#define COVARIANCE_MATRIX_H_

// ================================================================

#include <iostream>
#include <cassert>
#include <vector>
#include <algorithm>

// ================================================================

namespace zeno {
  
/// A covariance matrix.  Variables in the matrix are accessed by user-defined
/// ID numbers.
///
template <class T>
class CovarianceMatrix {
public:
  CovarianceMatrix(const CovarianceMatrix<T> & original);
  CovarianceMatrix();

  ~CovarianceMatrix();

  void add(unsigned int idToAdd, T variance);

  void copy(unsigned int sourceId, unsigned int destId);

  void clear(unsigned int idToClear);

  void remove(unsigned int idToRemove);

  void propagate(unsigned int toId, 
		 unsigned int fromIdA, T df_dA);

  void propagate(unsigned int toId,
		 unsigned int fromIdA, T df_dA,
		 unsigned int fromIdB, T df_dB);

  T getCovariance(unsigned int idA, 
		  unsigned int idB) const;

  void setCovariance(unsigned int idA, 
		     unsigned int idB, 
		     T covariance);

  unsigned int getSize() const;

  void print() const;

private:
  unsigned int idToMatrixIndex(unsigned int id) const;

  T getCovarianceFromIndexes(unsigned int indexA, 
			     unsigned int indexB) const;

  void setCovarianceFromIndexes(unsigned int indexA, 
				unsigned int indexB,
				T covariance);

  std::vector<unsigned int> activeIds;

  std::vector<std::vector<T> > matrix;
};

// ================================================================

/// Copy constructor.
///
template <class T>
CovarianceMatrix<T>::
CovarianceMatrix(const CovarianceMatrix<T> & original)
  : activeIds(original.activeIds),
    matrix(original.matrix)
{

}

/// Construct an empty covariance matrix.
///
template <class T>
CovarianceMatrix<T>::
CovarianceMatrix()
  : activeIds(),
    matrix()
{

}

template <class T>
CovarianceMatrix<T>::
~CovarianceMatrix()
{

}

/// Add a new variable with the given ID and variance.  Covariances of this
/// variable default to 0.
///
template <class T>
void
CovarianceMatrix<T>::
add(unsigned int idToAdd, T variance)
{
  activeIds.push_back(idToAdd);

  matrix.emplace_back(activeIds.size(), 0);

  matrix.back().back() = variance;
}

/// Copy the variance and covariances from one variable in the matrix to
/// another.  The covariance bewteen the variables and the variance of the
/// dest variable are set to the variance of the source variable.
///
template <class T>
void
CovarianceMatrix<T>::
copy(unsigned int sourceId, unsigned int destId) 
{
  unsigned int sourceIndex = idToMatrixIndex(sourceId);
  unsigned int destIndex   = idToMatrixIndex(destId);

  T sourceVariance = getCovarianceFromIndexes(sourceIndex, sourceIndex);

  setCovarianceFromIndexes(sourceIndex, destIndex, sourceVariance); 

  for (unsigned int i = 0; i < activeIds.size(); ++i) {
    T sourceCovariance = getCovarianceFromIndexes(i, sourceIndex);

    setCovarianceFromIndexes(i, destIndex, sourceCovariance);
  }
}

/// Set the variance and covariances of the given variable to zero.
///
template <class T>
void
CovarianceMatrix<T>::
clear(unsigned int idToClear)
{
  unsigned int indexToClear = idToMatrixIndex(idToClear);

  for (unsigned int i = 0; i < activeIds.size(); ++i) {
    setCovarianceFromIndexes(i, indexToClear, 0);
  }
}

/// Remove the given variable from the matrix.
///
template <class T>
void
CovarianceMatrix<T>::
remove(unsigned int idToRemove)
{
  unsigned int idIndex = idToMatrixIndex(idToRemove);

  for (unsigned int i = idIndex + 1; i < activeIds.size(); ++i) {
    matrix.at(i).erase(matrix.at(i).begin() + idIndex);
  }

  matrix.erase(matrix.begin() + idIndex);

  activeIds.erase(activeIds.begin() + idIndex);
}

template <class T>
T 
CovarianceMatrix<T>::
getCovariance(unsigned int idA, 
	      unsigned int idB) const
{
  unsigned int indexA = idToMatrixIndex(idA);
  unsigned int indexB = idToMatrixIndex(idB);

  return getCovarianceFromIndexes(indexA, indexB);
}

template <class T>
void 
CovarianceMatrix<T>::
setCovariance(unsigned int idA, 
	      unsigned int idB, 
	      T covariance)
{
  unsigned int indexA = idToMatrixIndex(idA);
  unsigned int indexB = idToMatrixIndex(idB);

  return setCovarianceFromIndexes(indexA, indexB, covariance);
}

template <class T>
unsigned int
CovarianceMatrix<T>::
getSize() const
{
  return activeIds.size();
}

template <class T>
void
CovarianceMatrix<T>::
print() const
{
  for (unsigned int i = 0; i < matrix.size(); ++ i) {
    for (unsigned int j = 0; j <= i; ++ j) {
      std::cout << matrix[i][j] << " ";
    }

    std::cout << std::endl;
  } 
}

/// Perform fist-order propagation of covariance from a variable A to a new
/// variable defined as f(A), given the derivative df/dA.
///
template <class T>
void 
CovarianceMatrix<T>::
propagate(unsigned int toId, 
	  unsigned int fromIdA, T df_dA)
{
  unsigned int indexA = idToMatrixIndex(fromIdA);

  unsigned int indexX = idToMatrixIndex(toId);

  for (unsigned int indexY = 0; indexY < matrix.size(); ++ indexY) {
    T covariance = 
      getCovarianceFromIndexes(indexY, indexA) * df_dA;

    setCovarianceFromIndexes(indexX, indexY, covariance);
  }
}

/// Perform fist-order propagation of covariance from variables A and B to a new
/// variable defined as f(A, B), given the partial derivatives df/dA and df/dB.
///
template <class T>
void 
CovarianceMatrix<T>::
propagate(unsigned int toId,
	  unsigned int fromIdA, T df_dA,
	  unsigned int fromIdB, T df_dB)
{
  unsigned int indexA = idToMatrixIndex(fromIdA);
  unsigned int indexB = idToMatrixIndex(fromIdB);

  unsigned int indexX = idToMatrixIndex(toId);

  for (unsigned int indexY = 0; indexY < matrix.size(); ++ indexY) {
    T covariance = 
      getCovarianceFromIndexes(indexY, indexA) * df_dA +
      getCovarianceFromIndexes(indexY, indexB) * df_dB;

    setCovarianceFromIndexes(indexX, indexY, covariance);
  }
}

/// Convert the given variable ID to an internal matrix index.
///
template <class T>
unsigned int 
CovarianceMatrix<T>::
idToMatrixIndex(unsigned int id) const
{
  std::vector<unsigned int>::const_iterator indexIt;

  indexIt = std::lower_bound(activeIds.begin(),
			     activeIds.end(),
			     id);

  unsigned int index = indexIt - activeIds.begin();

  return index;
}

template <class T>
T 
CovarianceMatrix<T>::
getCovarianceFromIndexes(unsigned int indexA, 
			 unsigned int indexB) const
{
  if (indexA < indexB) {
    std::swap(indexA, indexB);
  }

  T covariance = matrix.at(indexA).at(indexB);

  return covariance;
}

template <class T>
void 
CovarianceMatrix<T>::
setCovarianceFromIndexes(unsigned int indexA, 
			 unsigned int indexB,
			 T covariance)
{
  if (indexA < indexB) {
    std::swap(indexA, indexB);
  }

  matrix.at(indexA).at(indexB) = covariance;
}

}

#endif  // #ifndef COVARIANCE_MATRIX_H_

