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

#ifndef VECTOR3_H_
#define VECTOR3_H_

// ================================================================

#include <cstdio>
#include <cmath>
#include <cstring>
#include <cassert>
#include <ostream>
//#include <xmmintrin.h>

/// Represents a 3-dimensional vector.
///
template <class eleT>
class Vector3
{
 public:
  Vector3();
  Vector3(const Vector3<eleT> & b);
  Vector3(eleT newX, eleT newY, eleT newZ);

  eleT getX() const;
  eleT getY() const;
  eleT getZ() const;

  eleT getI() const;
  eleT getJ() const;
  eleT getK() const;

  void getXYZ(eleT & newX, eleT & newY, eleT & newZ) const;

  eleT get(int dimension) const;

  eleT getMin() const;
  eleT getMax() const;

  eleT dot(const Vector3<eleT> & b) const;
  Vector3<eleT> cross(const Vector3<eleT> & b) const;

  eleT getMagnitude() const;
  eleT getMagnitudeSqr() const;

  eleT getL1Magnitude() const;
  eleT getL1MagnitudeSqr() const;

  void setX(eleT newX);
  void setY(eleT newY);
  void setZ(eleT newZ);

  void setI(eleT newI);
  void setJ(eleT newJ);
  void setK(eleT newK);

  void setXYZ(eleT newX, eleT newY, eleT newZ);
  void setXYZ(eleT * newComponents);

  void set(int dimension, eleT value);

  void capBelow(const Vector3<eleT> & minValue);
  void capAbove(const Vector3<eleT> & maxValue);

  void normalize();

  Vector3<eleT> normalized() const;

  void translate(eleT deltaX, eleT deltaY, eleT deltaZ);

  void rotate(eleT cosTheta, eleT sinTheta, const Vector3<eleT> & axis);

  //in [col][row] format
  void leftMatrixMult3x3(eleT matrix[4][4]);

  void sort();

  Vector3<eleT> roundDown() const;
  Vector3<eleT> roundUp() const;

  void add(int dimension, eleT b);
  void mult(int dimension, eleT b);

  eleT & operator[](std::size_t idx);
  eleT const & operator[](std::size_t idx) const;

  template <class otherEleT> Vector3<eleT>
  operator+(const Vector3<otherEleT> & b) const;

  template <class otherEleT> Vector3<eleT>
  operator-(const Vector3<otherEleT> & b) const;

  template <class otherEleT> Vector3<eleT>
  operator*(const Vector3<otherEleT> & b) const;

  template <class otherEleT> Vector3<eleT>
  operator/(const Vector3<otherEleT> & b) const;

  Vector3<eleT> operator+(eleT b) const;
  Vector3<eleT> operator-(eleT b) const;
  Vector3<eleT> operator*(eleT b) const;
  Vector3<eleT> operator/(eleT b) const;
 
  Vector3<eleT> & operator+=(const Vector3<eleT> & rhs);
  Vector3<eleT> & operator-=(const Vector3<eleT> & rhs);
  Vector3<eleT> & operator*=(const Vector3<eleT> & rhs);
  Vector3<eleT> & operator/=(const Vector3<eleT> & rhs);

  Vector3<eleT> & operator+=(eleT b);
  Vector3<eleT> & operator-=(eleT b);
  Vector3<eleT> & operator*=(eleT b);
  Vector3<eleT> & operator/=(eleT b);

  template <class otherEleT>
  Vector3<eleT> & operator=(const Vector3<otherEleT> & rhs);

  template <class otherEleT>
  bool operator==(const Vector3<otherEleT> & rhs) const;

 private:
  void swap(eleT & a, eleT & b);
  
  eleT invSqrt(eleT x);

  eleT components[3];
};

template <class eleT>
std::ostream & operator<<(std::ostream & os, const Vector3<eleT> & rhs)
{
  os << "< " 
     << rhs.getX() << " , " << rhs.getY() << " , " << rhs.getZ() 
     << " >";

  return os;
}

template <class eleT>
Vector3<eleT>::Vector3()
{
  setXYZ(0, 0, 0);
}

template <class eleT>
Vector3<eleT>::Vector3(const Vector3<eleT> & b)
{
  setXYZ(b.getX(), b.getY(), b.getZ());
}

template <class eleT>
Vector3<eleT>::Vector3(eleT newX, eleT newY, eleT newZ)
{
  setXYZ(newX, newY, newZ);
}

template <class eleT>
eleT Vector3<eleT>::getX() const
{
  return components[0];
}

template <class eleT>
eleT Vector3<eleT>::getY() const
{
  return components[1];
}

template <class eleT>
eleT Vector3<eleT>::getZ() const
{
  return components[2];
}

template <class eleT>
eleT Vector3<eleT>::getI() const
{
  return getX();
}

template <class eleT>
eleT Vector3<eleT>::getJ() const
{
  return getY();
}

template <class eleT>
eleT Vector3<eleT>::getK() const
{
  return getZ();
}

template <class eleT>
void Vector3<eleT>::getXYZ(eleT & newX, eleT & newY, eleT & newZ) const
{
  newX = getX();
  newY = getY();
  newZ = getZ();
}

template <class eleT>
eleT Vector3<eleT>::get(int dimension) const
{
  return components[dimension];

  // switch (dimension) {
  //   case 0: return getX();
  //   case 1: return getY();
  //   case 2: return getZ();
  //   default:
  //     printf("Error: tried to return invalid dimension %i from vector\n",
  //            dimension);
  //     return 0;
  // }
}

template <class eleT>
eleT Vector3<eleT>::getMin() const
{
  if (get(0) < get(1)) {
    if (get(0) < get(2)) {
      return get(0);
    }
    else {
      return get(2);
    }
  }
  else {
    if (get(1) < get(2)) {
      return get(1);
    }
    else {
      return get(2);
    }
  }
}

template <class eleT>
eleT Vector3<eleT>::getMax() const
{
  if (get(0) > get(1)) {
    if (get(0) > get(2)) {
      return get(0);
    }
    else {
      return get(2);
    }
  }
  else {
    if (get(1) > get(2)) {
      return get(1);
    }
    else {
      return get(2);
    }
  }
}

template <class eleT>
eleT Vector3<eleT>::getMagnitude() const
{
  return sqrt(getX()*getX() + getY()*getY() + getZ()*getZ());
}

template <class eleT>
eleT Vector3<eleT>::getMagnitudeSqr() const
{
  return getX()*getX() + getY()*getY() + getZ()*getZ();
}

template <class eleT>
eleT Vector3<eleT>::getL1Magnitude() const
{
  return fabs(getX()) + fabs(getY()) + fabs(getZ());
}

template <class eleT>
eleT Vector3<eleT>::getL1MagnitudeSqr() const
{
  return pow(fabs(getX()) + fabs(getY()) + fabs(getZ()), 2);
}

template <class eleT>
void Vector3<eleT>::setX(eleT newX)
{
  components[0] = newX;
}

template <class eleT>
void Vector3<eleT>::setY(eleT newY)
{
  components[1] = newY;
}

template <class eleT>
void Vector3<eleT>::setZ(eleT newZ)
{
  components[2] = newZ;
}

template <class eleT>
void Vector3<eleT>::setI(eleT newI)
{
  setX(newI);
}

template <class eleT>
void Vector3<eleT>::setJ(eleT newJ)
{
  setY(newJ);
}

template <class eleT>
void Vector3<eleT>::setK(eleT newK)
{
  setZ(newK);
}

template <class eleT>
void Vector3<eleT>::setXYZ(eleT newX, eleT newY, eleT newZ)
{
  setX(newX);
  setY(newY);
  setZ(newZ);
}

template <class eleT>
void Vector3<eleT>::setXYZ(eleT *newComponents)
{
  memcpy(components, newComponents, 3*sizeof(eleT));
}

template <class eleT>
void Vector3<eleT>::set(int dimension, eleT value)
{
  components[dimension] = value;

  // switch (dimension) {
  //   case 0: setX(value); break;
  //   case 1: setY(value); break;
  //   case 2: setZ(value); break;
  //   default:
  //     printf("Error: tried to set invalid dimension %i in vector\n",
  //            dimension);
  //     break;
  // }
}

template <class eleT>
Vector3<eleT> Vector3<eleT>::roundDown() const
{
  return Vector3<eleT>(floor(getX()), floor(getY()), floor(getZ()));
}

template <class eleT>
Vector3<eleT> Vector3<eleT>::roundUp() const
{
  return Vector3<eleT>(ceil(getX()), ceil(getY()), ceil(getZ()));
}

template <class eleT>
void Vector3<eleT>::add(int dimension, eleT b)
{
  set(dimension, get(dimension) + b);
}

template <class eleT>
void Vector3<eleT>::mult(int dimension, eleT b)
{
  set(dimension, get(dimension) * b);
}

template <class eleT>
eleT &
Vector3<eleT>::operator[](std::size_t idx) {
  assert(idx >= 0 && idx < 3);

  return components[idx];
}

template <class eleT>
eleT const &
Vector3<eleT>::operator[](std::size_t idx) const {
  assert(idx >= 0 && idx < 3);

  return components[idx];
}

template <class eleT>
template <class otherEleT>
Vector3<eleT> Vector3<eleT>::operator+(const Vector3<otherEleT> & b) const
{
  Vector3<eleT> result(getX() + b.getX(), getY() + b.getY(),
                       getZ() + b.getZ());

  return result;
}

template <class eleT>
template <class otherEleT>
Vector3<eleT> Vector3<eleT>::operator-(const Vector3<otherEleT> & b) const
{
  Vector3<eleT> result(getX() - b.getX(), getY() - b.getY(),
                       getZ() - b.getZ());

  return result;
}

template <class eleT>
template <class otherEleT>
Vector3<eleT> Vector3<eleT>::operator*(const Vector3<otherEleT> & b) const
{
  Vector3<eleT> result(getX() * b.getX(), getY() * b.getY(),
                       getZ() * b.getZ());

  return result;
}

template <class eleT>
template <class otherEleT>
Vector3<eleT> Vector3<eleT>::operator/(const Vector3<otherEleT> & b) const
{
  Vector3<eleT> result(getX() / b.getX(), getY() / b.getY(),
                       getZ() / b.getZ());

  return result;
}

template <class eleT>
Vector3<eleT> Vector3<eleT>::operator+(eleT b) const
{
  Vector3<eleT> result(getX()+b, getY()+b, getZ()+b);

  return result;
}

template <class eleT>
Vector3<eleT> Vector3<eleT>::operator-(eleT b) const
{
  Vector3<eleT> result(getX()-b, getY()-b, getZ()-b);

  return result;
}

template <class eleT>
Vector3<eleT> Vector3<eleT>::operator*(eleT b) const
{
  Vector3<eleT> result(getX()*b, getY()*b, getZ()*b);

  return result;
}

template <class eleT>
Vector3<eleT> Vector3<eleT>::operator/(eleT b) const
{
  Vector3<eleT> result(getX()/b, getY()/b, getZ()/b);

  return result;
}

template <class eleT>
Vector3<eleT> operator+(eleT a, Vector3<eleT> b) 
{
  Vector3<eleT> result(a+b.getX(), a+b.getY(), a+b.getZ());

  return result;
}

template <class eleT>
Vector3<eleT> operator-(eleT a, Vector3<eleT> b) 
{
  Vector3<eleT> result(a-b.getX(), a-b.getY(), a-b.getZ());

  return result;
}

template <class eleT>
Vector3<eleT> operator*(eleT a, Vector3<eleT> b) 
{
  Vector3<eleT> result(a*b.getX(), a*b.getY(), a*b.getZ());

  return result;
}

template <class eleT>
Vector3<eleT> operator/(eleT a, Vector3<eleT> b) 
{
  Vector3<eleT> result(a/b.getX(), a/b.getY(), a/b.getZ());

  return result;
}

template <class eleT>
Vector3<eleT> & Vector3<eleT>::operator+=(const Vector3<eleT> & rhs)
{
  setXYZ(getX() + rhs.getX(), getY() + rhs.getY(), getZ() + rhs.getZ());

  return *this;
}

template <class eleT>
Vector3<eleT> & Vector3<eleT>::operator-=(const Vector3<eleT> & rhs)
{
  setXYZ(getX() - rhs.getX(), getY() - rhs.getY(), getZ() - rhs.getZ());

  return *this;
}

template <class eleT>
Vector3<eleT> & Vector3<eleT>::operator*=(const Vector3<eleT> & rhs)
{
  setXYZ(getX() * rhs.getX(), getY() * rhs.getY(), getZ() * rhs.getZ());

  return *this;
}

template <class eleT>
Vector3<eleT> & Vector3<eleT>::operator/=(const Vector3<eleT> & rhs)
{
  setXYZ(getX() / rhs.getX(), getY() / rhs.getY(), getZ() / rhs.getZ());

  return *this;
}

template <class eleT>
Vector3<eleT> & Vector3<eleT>::operator+=(eleT b)
{
  setXYZ(getX()+b, getY()+b, getZ()+b);

  return *this;
}

template <class eleT>
Vector3<eleT> & Vector3<eleT>::operator-=(eleT b)
{
  setXYZ(getX()-b, getY()-b, getZ()-b);

  return *this;
}

template <class eleT>
Vector3<eleT> & Vector3<eleT>::operator*=(eleT b)
{
  setXYZ(getX()*b, getY()*b, getZ()*b);

  return *this;
}

template <class eleT>
Vector3<eleT> & Vector3<eleT>::operator/=(eleT b)
{
  setXYZ(getX()/b, getY()/b, getZ()/b);

  return *this;
}

template <class eleT>
template <class otherEleT>
Vector3<eleT> & Vector3<eleT>::operator=(const Vector3<otherEleT> & rhs)
{
  setXYZ((eleT)rhs.getX(), (eleT)rhs.getY(), (eleT)rhs.getZ());

  return *this;
}

template <class eleT>
template <class otherEleT>
bool Vector3<eleT>::operator==(const Vector3<otherEleT> & rhs) const
{
  return ((getX() == rhs.getX()) &&
	  (getY() == rhs.getY()) &&
	  (getZ() == rhs.getZ()));
}

template <class eleT>
void Vector3<eleT>::translate(eleT deltaX, eleT deltaY, eleT deltaZ)
{
  setXYZ(getX() + deltaX, getY() + deltaY, getZ() + deltaZ);
}

template <class eleT>
void Vector3<eleT>::rotate(eleT cosTheta, eleT sinTheta,
			   const Vector3<eleT> & axis)
{
  Vector3<eleT> u(axis);
  u.normalize();

  eleT matrix[4][4];

  matrix[0][0] = cosTheta + pow(u.getX(), 2)*(1 - cosTheta); 
  matrix[0][1] = u.getX()*u.getY()*(1 - cosTheta) - u.getZ()*sinTheta;
  matrix[0][2] = u.getX()*u.getZ()*(1 - cosTheta) + u.getY()*sinTheta;
  matrix[0][3] = 0;

  matrix[1][0] = u.getY()*u.getX()*(1 - cosTheta) + u.getZ()*sinTheta;
  matrix[1][1] = cosTheta + pow(u.getY(), 2)*(1 - cosTheta);
  matrix[1][2] = u.getY()*u.getZ()*(1 - cosTheta) - u.getX()*sinTheta;
  matrix[1][3] = 0;

  matrix[2][0] = u.getZ()*u.getX()*(1 - cosTheta) - u.getY()*sinTheta;
  matrix[2][1] = u.getZ()*u.getY()*(1 - cosTheta) + u.getX()*sinTheta;
  matrix[2][2] = cosTheta + pow(u.getZ(), 2)*(1 - cosTheta);
  matrix[2][3] = 0;

  matrix[3][0] = 0;
  matrix[3][1] = 0;
  matrix[3][2] = 0;
  matrix[3][3] = 1;

  leftMatrixMult3x3(matrix);
}

template <class eleT>
eleT Vector3<eleT>::dot(const Vector3<eleT> & b) const
{
  return (getX()*b.getX() + getY()*b.getY() + getZ()*b.getZ());
}

template <class eleT>
Vector3<eleT> Vector3<eleT>::cross(const Vector3<eleT> & b) const
{
  eleT crossX = getY()*b.getZ() - getZ()*b.getY();
  eleT crossY = getZ()*b.getX() - getX()*b.getZ();
  eleT crossZ = getX()*b.getY() - getY()*b.getX();

  Vector3<eleT> result(crossX, crossY, crossZ);

  return result;
}

template <class eleT>
void Vector3<eleT>::capBelow(const Vector3<eleT> & minValue)
{
  if (getX() < minValue.getX()) setX(minValue.getX());
  if (getY() < minValue.getY()) setY(minValue.getY());
  if (getZ() < minValue.getZ()) setZ(minValue.getZ());
}

template <class eleT>
void Vector3<eleT>::capAbove(const Vector3<eleT> & maxValue)
{
  if (getX() > maxValue.getX()) setX(maxValue.getX());
  if (getY() > maxValue.getY()) setY(maxValue.getY());
  if (getZ() > maxValue.getZ()) setZ(maxValue.getZ());
}

template <class eleT>
void Vector3<eleT>::normalize()
{
  eleT magnitudeSqr = getMagnitudeSqr();

  if (magnitudeSqr == 0)
    return;

  eleT invMagnitude = invSqrt(magnitudeSqr);

  setXYZ(getX()*invMagnitude, getY()*invMagnitude, getZ()*invMagnitude);
}

template <class eleT>
Vector3<eleT> Vector3<eleT>::normalized() const
{
  Vector3<eleT> normalized(*this);
  normalized.normalize();

  return normalized;
}

//in [col][row] format
template <class eleT>
void Vector3<eleT>::leftMatrixMult3x3(eleT matrix[4][4])
{
  eleT _x = getX();
  eleT _y = getY();
  eleT _z = getZ();

  setX(_x * matrix[0][0] + _y * matrix[1][0] + _z * matrix[2][0]);
  setY(_x * matrix[0][1] + _y * matrix[1][1] + _z * matrix[2][1]);
  setZ(_x * matrix[0][2] + _y * matrix[1][2] + _z * matrix[2][2]);
  /*
    setX(_x * matrix[0][0] + _y * matrix[0][1] + _z * matrix[0][2]);
    setY(_x * matrix[1][0] + _y * matrix[1][1] + _z * matrix[1][2]);
    setZ(_x * matrix[2][0] + _y * matrix[2][1] + _z * matrix[2][2]);
  */
}

template <class eleT>
void Vector3<eleT>::sort()
{
  if(components[1]<components[0]) 
    swap(components[0],components[1]);
  
  if(components[2]<components[0]) 
    swap(components[0],components[2]);
  
  if(components[2]<components[1]) 
    swap(components[1],components[2]);
}

template <class eleT>
void Vector3<eleT>::swap(eleT & a, eleT & b)
{
  eleT temp = a;
  a = b;
  b = temp;
}

template <class eleT>
eleT Vector3<eleT>::invSqrt(eleT x)
{
  return 1. / sqrt(x);

  //return _mm_cvtss_si32( _mm_rsqrt_ss( _mm_set_ss( x ) ) );

  //union {
  //	float xFloat;
  //	int xInt;
  //};
  //xFloat = x;
  //xInt = 0x5f3759df - (xInt >> 1); //compute initial guess for Newton iteration
  //x = xFloat;
  //float xhalf = 0.5f*x;
  //x = x * (1.5f - xhalf * x * x); //Newton iteration step
  //x = x * (1.5f - xhalf * x * x);
  //x = x * (1.5f - xhalf * x * x);
  //return x;
}

// ================================================================

#endif  // #ifndef VECTOR3_H_

