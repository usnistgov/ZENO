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
// Time-stamp: <2016-09-19 15:24:35 dcj>
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
  void cubicsolver(eleT coeff[4], 
		   eleT roots[3]) const;

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
  eleT coeffs[4];
  eleT roots[3];

  assert(get(1, 0) == get(0, 1));
  assert(get(2, 0) == get(0, 2));
  assert(get(2, 1) == get(1, 2));

  // compute Rg tensor
  eleT Sxx=get(0, 0);
  eleT Syy=get(1, 1);
  eleT Szz=get(2, 2);
  eleT Sxy=get(0, 1);
  eleT Sxz=get(0, 2);
  eleT Syz=get(1, 2);

  // Compute coefficents
  coeffs[0]=1.;
  coeffs[1]=-1.*(Sxx+Syy+Szz);
  coeffs[2]=Sxx*Syy+Sxx*Szz+Syy*Szz-Sxy*Sxy-Sxz*Sxz-Syz*Syz;
  coeffs[3]=Sxz*Sxz*Syy+Sxy*Sxy*Szz+Syz*Syz*Sxx-2.*Sxy*Sxz*Syz-Sxx*Syy*Szz;

  // Find eigenvalues
  cubicsolver(coeffs,roots);

  eigenValues.setXYZ(roots[0], roots[1], roots[2]);
}

template <class eleT>
void Matrix3x3<eleT>::cubicsolver(eleT coeff[4], 
				  eleT roots[3]) const
{
// solve the cubic equation and sort roots

  eleT delta0 = coeff[1]*coeff[1]-3.*coeff[0]*coeff[2];
  eleT delta1 = 2.*pow(coeff[1],3.)-9.*coeff[0]*coeff[1]*coeff[2]+27.*coeff[0]*coeff[0]*coeff[3];

  eleT ththeta = pow(delta0,3.)-pow(delta1,2.)/4.;
  if(ththeta < 1e-15) {ththeta=0;}
  ththeta = atan2(sqrt(ththeta),delta1/2.)/3.;

  eleT cosv = cos(ththeta);
  eleT sinv = sin(ththeta);

  roots[0] = -1./(3.*coeff[0])*(coeff[1]+2.*sqrt(delta0)*cosv);
  roots[1] = -1./(3.*coeff[0])*(coeff[1]+sqrt(delta0)*(-1.*cosv-sqrt(3.)*sinv));
  roots[2] = -1./(3.*coeff[0])*(coeff[1]+sqrt(delta0)*(-1.*cosv+sqrt(3.)*sinv));

  // sort the roots
  if(roots[1]<roots[0]) 
    std::swap(roots[0],roots[1]);
  if(roots[2]<roots[0]) 
    std::swap(roots[0],roots[2]);
  if(roots[2]<roots[1]) 
    std::swap(roots[1],roots[2]);
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
