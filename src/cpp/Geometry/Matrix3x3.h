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
// Date:    Tue Apr 16 12:20:59 2013 EDT
//
// Time-stamp: <2017-03-24 18:05:28 dcj>
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

  const int maxJacobiIter = 1000;

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
	
	if ((fabs(workEigenVals[row]) + fabs(workMatrix[row][col]) ==
	     fabs(workEigenVals[row])) &&
	    
	    (fabs(workEigenVals[col]) + fabs(workMatrix[row][col]) ==
	     fabs(workEigenVals[col]))) {
	  
	  workMatrix[row][col] = 0.;
	}
	else {  
	  // Calculate Jacobi transformation

	  eleT tan_phi = eleT();
	  
	  eleT h = workEigenVals[col] - workEigenVals[row];
	  
	  if (fabs(h) + fabs(workMatrix[row][col]) ==
	      fabs(h)) {
	    
	    tan_phi = workMatrix[row][col] / h;
	  }
	  else {
	    eleT phi = 0.5 * h / workMatrix[row][col];
	    
	    if (phi < 0.)
	      tan_phi = -1. / (sqrt(1. + (phi*phi)) - phi);
	    else
	      tan_phi = 1. / (sqrt(1. + (phi*phi)) + phi);
	  }
	    
	  eleT cos_phi = 1. / sqrt(1. + (tan_phi*tan_phi));
	  eleT sin_phi = tan_phi * cos_phi;

	  // Apply Jacobi transformation

	  workEigenVals[row] -= tan_phi * workMatrix[row][col];
	  workEigenVals[col] += tan_phi * workMatrix[row][col];
	  
	  workMatrix[row][col] = 0.;
	    
	  for (int i = 0; i < row; i++) {
	    tan_phi = workMatrix[i][row];
	    workMatrix[i][row] = cos_phi*tan_phi - sin_phi*workMatrix[i][col];
	    workMatrix[i][col] = sin_phi*tan_phi + cos_phi*workMatrix[i][col];
	  }
	    
	  for (int i = row + 1; i < col; i++) {
	    tan_phi = workMatrix[row][i];
	    workMatrix[row][i] = cos_phi*tan_phi - sin_phi*workMatrix[i][col];
	    workMatrix[i][col] = sin_phi*tan_phi + cos_phi*workMatrix[i][col];
	  }
	    
	  for (int i = col + 1; i < dim; i++) {
	    tan_phi = workMatrix[row][i];
	    workMatrix[row][i] = cos_phi*tan_phi - sin_phi*workMatrix[col][i];
	    workMatrix[col][i] = sin_phi*tan_phi + cos_phi*workMatrix[col][i];
	  }
	}
      }
    }

    ++ numJacobiIter;

    assert(numJacobiIter < maxJacobiIter);
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

// ================================================================

// Local Variables:
// time-stamp-line-limit: 30
// mode: c++
// End:
