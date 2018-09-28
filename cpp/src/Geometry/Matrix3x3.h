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
// Created: Tue Apr 16 12:20:59 2013 EDT
//
// ================================================================

#ifndef MATRIX3X3_H_
#define MATRIX3X3_H_

// ================================================================

#include <cstdio>
#include <cmath>
#include <cstring>
#include <cassert>
#include <ostream>
#include <algorithm>

#include "Vector3.h"

/// Represents a 3x3 matrix.
///
template <class eleT>
class Matrix3x3
{
 public:
  Matrix3x3();
  Matrix3x3(const Matrix3x3<eleT> & matrix);
  Matrix3x3(const eleT * matrix);
  Matrix3x3(eleT matrix00, eleT matrix01, eleT matrix02,
	    eleT matrix10, eleT matrix11, eleT matrix12,
	    eleT matrix20, eleT matrix21, eleT matrix22);

  eleT get(int row, int col) const;
  eleT get(int component) const;

  void set(int row, int col, eleT value);
  void set(int component, eleT value);

  void symmetrize();

  void getEigenValues(Vector3<eleT> & eigenValues) const;

  void addRow(int row, const Vector3<eleT> & b);
  void addCol(int col, const Vector3<eleT> & b);

  void multRow(int row, const Vector3<eleT> & b);
  void multCol(int col, const Vector3<eleT> & b);

  template <class otherEleT> Matrix3x3<eleT>
  operator+(const Matrix3x3<otherEleT> & b) const;

  template <class otherEleT> Matrix3x3<eleT>
  operator-(const Matrix3x3<otherEleT> & b) const;

  template <class otherEleT> Matrix3x3<eleT>
  operator*(const Matrix3x3<otherEleT> & b) const;

  template <class otherEleT> Matrix3x3<eleT>
  operator/(const Matrix3x3<otherEleT> & b) const;

  Matrix3x3<eleT> operator+(eleT b) const;
  Matrix3x3<eleT> operator-(eleT b) const;
  Matrix3x3<eleT> operator*(eleT b) const;
  Matrix3x3<eleT> operator/(eleT b) const;

  Matrix3x3<eleT> & operator+=(const Matrix3x3<eleT> & rhs);
  Matrix3x3<eleT> & operator-=(const Matrix3x3<eleT> & rhs);
  Matrix3x3<eleT> & operator*=(const Matrix3x3<eleT> & rhs);
  Matrix3x3<eleT> & operator/=(const Matrix3x3<eleT> & rhs);

  Matrix3x3<eleT> & operator+=(eleT b);
  Matrix3x3<eleT> & operator-=(eleT b);
  Matrix3x3<eleT> & operator*=(eleT b);
  Matrix3x3<eleT> & operator/=(eleT b);

  template <class otherEleT>
  Matrix3x3<eleT> & operator=(const Matrix3x3<otherEleT> & rhs);

 private:
  eleT components[3*3];
};

template <class eleT>
std::ostream & operator<<(std::ostream & os, const Matrix3x3<eleT> & rhs)
{
  os << "[ " 
     << rhs.get(0, 0) << " , " 
     << rhs.get(0, 1) << " , " 
     << rhs.get(0, 2) << " ,"
     << std::endl
     << "  "
     << rhs.get(1, 0) << " , " 
     << rhs.get(1, 1) << " , " 
     << rhs.get(1, 2) << " ,"
     << std::endl
     << "  "
     << rhs.get(2, 0) << " , " 
     << rhs.get(2, 1) << " , " 
     << rhs.get(2, 2) << " ]"
     << std::endl;

  return os;
}

template <class eleT>
Matrix3x3<eleT>::Matrix3x3()
{
  for (int i = 0; i < 3*3; i++) {
    components[i] = 0;
  }
}

template <class eleT>
Matrix3x3<eleT>::Matrix3x3(const Matrix3x3<eleT> & matrix)
{
  for (int i = 0; i < 3*3; i++) {
    components[i] = matrix.components[i];
  }
}

template <class eleT>
Matrix3x3<eleT>::Matrix3x3(const eleT * matrix)
{
  for (int i = 0; i < 3*3; i++) {
    components[i] = matrix[i];
  }
}

template <class eleT>
Matrix3x3<eleT>::Matrix3x3(eleT matrix00, eleT matrix01, eleT matrix02,
			   eleT matrix10, eleT matrix11, eleT matrix12,
			   eleT matrix20, eleT matrix21, eleT matrix22)
{
  set(0, 0, matrix00);
  set(0, 1, matrix01);
  set(0, 2, matrix02);
  set(1, 0, matrix10);
  set(1, 1, matrix11);
  set(1, 2, matrix12);
  set(2, 0, matrix20);
  set(2, 1, matrix21);
  set(2, 2, matrix22);
}

template <class eleT>
eleT Matrix3x3<eleT>::get(int row, int col) const
{
  assert(row >= 0 && row < 3);
  assert(col >= 0 && col < 3);

  return components[row*3 + col];
}

template <class eleT>
eleT Matrix3x3<eleT>::get(int component) const
{
  assert(component >= 0 && component < 9);

  return components[component];
}

template <class eleT>
void Matrix3x3<eleT>::set(int row, int col, eleT value)
{
  assert(row >= 0 && row < 3);
  assert(col >= 0 && col < 3);

  components[row*3 + col] = value;
}

template <class eleT>
void Matrix3x3<eleT>::set(int component, eleT value)
{
  assert(component >= 0 && component < 9);

  components[component] = value;
}

template <class eleT>
void Matrix3x3<eleT>::symmetrize()
{
  for (int row = 0; row < 3; row++) {
    for (int col = row + 1; col < 3; col++) {
      eleT average = (get(row, col) + get(col, row)) / (eleT)2;

      set(row, col, average);
      set(col, row, average);
    }
  }
}

template <class eleT>
void Matrix3x3<eleT>::getEigenValues(Vector3<eleT> & eigenValues) const
{
  // Jacobi
  
  // Joachim Kopp
  // Efficient numerical diagonalization of hermitian 3x3 matrices
  // Int. J. Mod. Phys. C 19 (2008) 523-548
  // arXiv.org: physics/0610206
  
  assert(get(1, 0) == get(0, 1));
  assert(get(2, 0) == get(0, 2));
  assert(get(2, 1) == get(1, 2));

  const int dim = 3;

  const int maxJacobiIter = 100;

  eleT workEigenVals[3];

  eleT workMatrix[3][3];

  for (int i = 0; i < dim; i++)
    for (int j = 0; j < dim; j++)
      workMatrix[i][j] = get(i, j);
  
  for (int i = 0; i < dim; i++)
    workEigenVals[i] = workMatrix[i][i];

  int numJacobiIter = 0;

  // Do Jacobi iterations
  while ((workMatrix[0][1] != 0.) ||
         (workMatrix[0][2] != 0.) ||
	 (workMatrix[1][2] != 0.)) {

    // Do sweep
    for (int row = 0; row < dim; row++) {
      for (int col = row + 1; col < dim; col++) {

	// Calculate Jacobi transformation

	eleT tanPhi = eleT();
	  
	eleT thetaNumer = workEigenVals[col] - workEigenVals[row];
	eleT thetaDenom = 2. * workMatrix[row][col];

	// if thetaNumer/thetaDenom might overflow
	if (fabs(thetaNumer) + fabs(thetaDenom) ==
	    fabs(thetaNumer)) {
	  
	  // tan(phi) = 1/(2*theta)
	  tanPhi = thetaDenom / (2. * thetaNumer);
	}
	else {
	  eleT theta = thetaNumer / thetaDenom;
	    
	  if (theta < 0.) {
	    tanPhi = -1. / (sqrt(1. + pow(theta, 2.)) - theta);
	  }
	  else {
	    tanPhi = 1. / (sqrt(1. + pow(theta, 2.)) + theta);
	  }
	}

	eleT cosPhi = 1. / sqrt(1. + pow(tanPhi, 2.));
	eleT sinPhi = tanPhi * cosPhi;

	// Apply Jacobi transformation

	workEigenVals[row] -= tanPhi * workMatrix[row][col];
	workEigenVals[col] += tanPhi * workMatrix[row][col];
	  
	workMatrix[row][col] = 0.;
	    
	for (int i = 0; i < row; i++) {
	  tanPhi = workMatrix[i][row];
	  workMatrix[i][row] = cosPhi * tanPhi - sinPhi * workMatrix[i][col];
	  workMatrix[i][col] = sinPhi * tanPhi + cosPhi * workMatrix[i][col];
	}
	    
	for (int i = row + 1; i < col; i++) {
	  tanPhi = workMatrix[row][i];
	  workMatrix[row][i] = cosPhi * tanPhi - sinPhi * workMatrix[i][col];
	  workMatrix[i][col] = sinPhi * tanPhi + cosPhi * workMatrix[i][col];
	}
	    
	for (int i = col + 1; i < dim; i++) {
	  tanPhi = workMatrix[row][i];
	  workMatrix[row][i] = cosPhi * tanPhi - sinPhi * workMatrix[col][i];
	  workMatrix[col][i] = sinPhi * tanPhi + cosPhi * workMatrix[col][i];
	}
      }
    }

    ++ numJacobiIter;

    if (numJacobiIter > maxJacobiIter) {
      std::cerr << std::endl
                << "*** Warning ***" << std::endl
		<< "Jacobi iterations did not converge in getEigenValues.  "
	        << "Eigen values may not be correct."
	        << std::endl;

      break;
    }
  }

  eigenValues.setXYZ(workEigenVals[0], workEigenVals[1], workEigenVals[2]);
  eigenValues.sort();
}

template <class eleT>
void Matrix3x3<eleT>::addRow(int row, const Vector3<eleT> & b)
{
  assert(row >= 0 && row < 3);

  for (int col = 0; col < 3; col++) {
    eleT sum = get(row, col) + b.get(col);

    set(row, col, sum);
  }
}

template <class eleT>
void Matrix3x3<eleT>::addCol(int col, const Vector3<eleT> & b)
{
  assert(col >= 0 && col < 3);

  for (int row = 0; row < 3; row++) {
    eleT sum = get(row, col) + b.get(row);

    set(row, col, sum);
  }
}

template <class eleT>
void Matrix3x3<eleT>::multRow(int row, const Vector3<eleT> & b)
{
  assert(row >= 0 && row < 3);

  for (int col = 0; col < 3; col++) {
    eleT product = get(row, col) * b.get(col);

    set(row, col, product);
  }
}

template <class eleT>
void Matrix3x3<eleT>::multCol(int col, const Vector3<eleT> & b)
{
  assert(col >= 0 && col < 3);

  for (int row = 0; row < 3; row++) {
    eleT product = get(row, col) * b.get(row);

    set(row, col, product);
  }
}

template <class eleT>
template <class otherEleT>
Matrix3x3<eleT> Matrix3x3<eleT>::operator+(const Matrix3x3<otherEleT> & b) const
{
  eleT result[3*3];

  for (int i = 0; i < 3*3; i++) {
    result[i] = components[i] + b.components[i];
  }

  return Matrix3x3<eleT>(result);
}

template <class eleT>
template <class otherEleT>
Matrix3x3<eleT> Matrix3x3<eleT>::operator-(const Matrix3x3<otherEleT> & b) const
{
  eleT result[3*3];

  for (int i = 0; i < 3*3; i++) {
    result[i] = components[i] - b.components[i];
  }

  return Matrix3x3<eleT>(result);
}

template <class eleT>
template <class otherEleT>
Matrix3x3<eleT> Matrix3x3<eleT>::operator*(const Matrix3x3<otherEleT> & b) const
{
  eleT result[3*3];

  for (int i = 0; i < 3*3; i++) {
    result[i] = components[i] * b.components[i];
  }

  return Matrix3x3<eleT>(result);
}

template <class eleT>
template <class otherEleT>
Matrix3x3<eleT> Matrix3x3<eleT>::operator/(const Matrix3x3<otherEleT> & b) const
{
  eleT result[3*3];

  for (int i = 0; i < 3*3; i++) {
    result[i] = components[i] / b.components[i];
  }

  return Matrix3x3<eleT>(result);
}

template <class eleT>
Matrix3x3<eleT> Matrix3x3<eleT>::operator+(eleT b) const
{
  eleT result[3*3];

  for (int i = 0; i < 3*3; i++) {
    result[i] = components[i] + b;
  }

  return Matrix3x3<eleT>(result);
}

template <class eleT>
Matrix3x3<eleT> Matrix3x3<eleT>::operator-(eleT b) const
{
  eleT result[3*3];

  for (int i = 0; i < 3*3; i++) {
    result[i] = components[i] - b;
  }

  return Matrix3x3<eleT>(result);
}

template <class eleT>
Matrix3x3<eleT> Matrix3x3<eleT>::operator*(eleT b) const
{
  eleT result[3*3];

  for (int i = 0; i < 3*3; i++) {
    result[i] = components[i] * b;
  }

  return Matrix3x3<eleT>(result);
}

template <class eleT>
Matrix3x3<eleT> Matrix3x3<eleT>::operator/(eleT b) const
{
  eleT result[3*3];

  for (int i = 0; i < 3*3; i++) {
    result[i] = components[i] / b;
  }

  return Matrix3x3<eleT>(result);
}

template <class eleT>
Matrix3x3<eleT> & Matrix3x3<eleT>::operator+=(const Matrix3x3<eleT> & rhs)
{
  for (int i = 0; i < 3*3; i++) {
    components[i] += rhs.components[i];
  }

  return *this;
}

template <class eleT>
Matrix3x3<eleT> & Matrix3x3<eleT>::operator-=(const Matrix3x3<eleT> & rhs)
{
  for (int i = 0; i < 3*3; i++) {
    components[i] -= rhs.components[i];
  }

  return *this;
}

template <class eleT>
Matrix3x3<eleT> & Matrix3x3<eleT>::operator*=(const Matrix3x3<eleT> & rhs)
{
  for (int i = 0; i < 3*3; i++) {
    components[i] *= rhs.components[i];
  }

  return *this;
}

template <class eleT>
Matrix3x3<eleT> & Matrix3x3<eleT>::operator/=(const Matrix3x3<eleT> & rhs)
{
  for (int i = 0; i < 3*3; i++) {
    components[i] /= rhs.components[i];
  }

  return *this;
}

template <class eleT>
Matrix3x3<eleT> & Matrix3x3<eleT>::operator+=(eleT b)
{
  for (int i = 0; i < 3*3; i++) {
    components[i] += b;
  }

  return *this;
}

template <class eleT>
Matrix3x3<eleT> & Matrix3x3<eleT>::operator-=(eleT b)
{
  for (int i = 0; i < 3*3; i++) {
    components[i] -= b;
  }

  return *this;
}

template <class eleT>
Matrix3x3<eleT> & Matrix3x3<eleT>::operator*=(eleT b)
{
  for (int i = 0; i < 3*3; i++) {
    components[i] *= b;
  }

  return *this;
}

template <class eleT>
Matrix3x3<eleT> & Matrix3x3<eleT>::operator/=(eleT b)
{
  for (int i = 0; i < 3*3; i++) {
    components[i] /= b;
  }

  return *this;
}

template <class eleT>
template <class otherEleT>
Matrix3x3<eleT> & Matrix3x3<eleT>::operator=(const Matrix3x3<otherEleT> & rhs)
{
  for (int i = 0; i < 3*3; i++) {
    otherEleT rhsComponent = rhs.get(i);
    set(i, (eleT)rhsComponent);
  }

  return *this;
}

// ================================================================

#endif  // #ifndef MATRIX3X3_H_

