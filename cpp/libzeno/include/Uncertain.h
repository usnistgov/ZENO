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
// Created: Mon Nov 16 15:04:44 2015 EDT
//
// ================================================================

#ifndef UNCERTAIN_H_
#define UNCERTAIN_H_

// ================================================================

#include <iostream>
#include <cmath>
#include <cassert>

#include "CovarianceMatrix.h"

#include "Geometry/Vector3.h"
#include "Geometry/Matrix3x3.h"

// ================================================================

namespace zeno {

template <class T>
class Uncertain;

template <class T>
Uncertain<T> fabs(const Uncertain<T> & A); 

template <class T>
Uncertain<T> log(const Uncertain<T> & A); 

template <class T>
Uncertain<T> exp(const Uncertain<T> & A);

template <class T>
Uncertain<T> pow(const Uncertain<T> & A, const Uncertain<T> & B);

template <class T>
Uncertain<T> sqrt(const Uncertain<T> & A);

template <class T>
Uncertain<T> sin(const Uncertain<T> & A);

template <class T>
Uncertain<T> cos(const Uncertain<T> & A);

template <class T>
Uncertain<T> atan2(const Uncertain<T> & A, const Uncertain<T> & B);

// ================================================================

/// Represents a value with uncertainty as a mean and variance.
/// Handles uncertainty propagation through arithmatic operations and
/// tracks covariance between variables.
///
template <class T>
class Uncertain {
public:
  Uncertain(const Uncertain<T> & original);
  Uncertain(T mean = 0, T variance = 0);

  ~Uncertain();

  T getMean() const;
  T getVariance() const;

  T getStdDev() const;

  Uncertain<T> & operator=(const Uncertain<T> & rhs);
  Uncertain<T> & operator=(const T & rhs);

  Uncertain<T> & operator+=(const Uncertain<T> & rhs);
  Uncertain<T> & operator+=(const T & rhs);

  Uncertain<T> & operator-=(const Uncertain<T> & rhs);
  Uncertain<T> & operator-=(const T & rhs);

  Uncertain<T> & operator*=(const Uncertain<T> & rhs);
  Uncertain<T> & operator*=(const T & rhs);

  Uncertain<T> & operator/=(const Uncertain<T> & rhs);
  Uncertain<T> & operator/=(const T & rhs);

  friend Uncertain fabs<T>(const Uncertain & A); 
  friend Uncertain log<T>(const Uncertain & A); 
  friend Uncertain exp<T>(const Uncertain & A);
  friend Uncertain pow<T>(const Uncertain & A, const Uncertain & B);
  friend Uncertain sqrt<T>(const Uncertain & A);
  friend Uncertain sin<T>(const Uncertain & A);
  friend Uncertain cos<T>(const Uncertain & A);
  friend Uncertain atan2<T>(const Uncertain & A, const Uncertain & B);

  static Vector3<Uncertain<T> > zip(Vector3<T> const & mean,
				    Vector3<T> const & variance);

  static Matrix3x3<Uncertain<T> > zip(Matrix3x3<T> const & mean,
				      Matrix3x3<T> const & variance);

  static void setCovariance(Uncertain<T> const & A, Uncertain<T> const & B, double);

private:
  static T firstOrderVariance(T varA, T df_dA);

  static T firstOrderVariance(T varA, T df_dA,
			      T varB, T df_dB);

  static T firstOrderVariance(T varA, T df_dA, 
			      T covAB,
			      T varB, T df_dB);

  static T secondOrderMean(T firstOrderMean, 
			   T varA, T d2f_dA2);

  static T secondOrderMean(T firstOrderMean, 
			   T varA, T d2f_dA2,
			   T varB, T d2f_dB2);

  static T secondOrderMean(T firstOrderMean, 
			   T varA, T d2f_dA2,
			   T covAB, T d2f_dAdB, 
			   T varB, T d2f_dB2);

  static T secondOrderVariance(T firstOrderVariance, 
			       T varA, T d2f_dA2);

  static T secondOrderVariance(T firstOrderVariance, 
			       T varA, T d2f_dA2,
			       T d2f_dAdB,  
			       T varB, T d2f_dB2);

  static T secondOrderVariance(T firstOrderVariance, 
			       T varA, T d2f_dA2,
			       T covAB, T d2f_dAdB, 
			       T varB, T d2f_dB2);

  void updateFirstOrder(unsigned int idA, T varA, T df_dA);

  void updateFirstOrder(unsigned int idA, T varA, T df_dA,
			unsigned int idB, T varB, T df_dB);

  void updateSecondOrder(T varA, T d2f_dA2);

  void updateSecondOrder(T varA, T d2f_dA2,
			 T d2f_dAdB,
			 T varB, T d2f_dB2);

  static const bool propagateUncertainty = true;
  static const bool useCovariance        = true;
  static const bool useSecondOrder       = false;

  static const bool debug = false;

  T mean;
  T variance;

  unsigned int id;

  static unsigned int nextId;

  static CovarianceMatrix<T> covarianceMatrix;
};

// ================================================================

template <class T>
unsigned int Uncertain<T>::nextId = 0;

template <class T>
CovarianceMatrix<T> Uncertain<T>::covarianceMatrix;

// ================================================================

/// Copy constructor.
///
template <class T>
Uncertain<T>::
Uncertain(const Uncertain<T> & original) 
  : mean(original.mean),
    variance(original.variance),
    id()
{
  if (useCovariance) {
    id = nextId ++;

    covarianceMatrix.add(id, 0);
    covarianceMatrix.copy(original.id, id);
  }

  if (debug) {
    std::cout << "Construct new " << id 
	      << " (" << covarianceMatrix.getSize() << ")" << std::endl;

    covarianceMatrix.print();
  }
}

/// Constructs a new Uncertain value with the given mean and variance.
///
template <class T>
Uncertain<T>::
Uncertain(T mean, T variance)
  : mean(mean),
    variance(variance),
    id() 
{
  if (useCovariance) {
    id = nextId ++;

    covarianceMatrix.add(id, variance);
  }

  if (debug) {
    std::cout << "Construct copy " << id
	      << " (" << covarianceMatrix.getSize() << ")" << std::endl;

    covarianceMatrix.print();
  }
}

template <class T>
Uncertain<T>::
~Uncertain() 
{
  if (useCovariance) {
    covarianceMatrix.remove(id);
  }

  if (debug) {
    std::cout << "Destruct " << id
	      << " (" << covarianceMatrix.getSize() << ")" << std::endl;

    covarianceMatrix.print();
  }
}

template <class T>
T
Uncertain<T>::
getMean() const 
{
  return mean;
}


template <class T>
T
Uncertain<T>::
getVariance() const 
{
  return variance;
}

template <class T>
T
Uncertain<T>::
getStdDev() const 
{
  return std::sqrt(variance);
}

template <class T>
Uncertain<T> & 
Uncertain<T>::
operator=(const Uncertain<T> & rhs)
{
  mean     = rhs.mean;
  variance = rhs.variance;

  if (useCovariance) {
    covarianceMatrix.copy(rhs.id, id);
  }

  if (debug) {
    std::cout << id << " = " << rhs.id << std::endl;

    covarianceMatrix.print();
  }

  return *this;
}

template <class T>
Uncertain<T> & 
Uncertain<T>::
operator=(const T & rhs)
{
  mean     = rhs;
  variance = 0;

  if (useCovariance) {
    covarianceMatrix.clear(id);
  }

  return *this;
}

template <class T>
Uncertain<T> & 
Uncertain<T>::
operator+=(const Uncertain<T> & rhs)
{
  const T meanA = this->mean;
  const T meanB = rhs.mean;

  const T varA = this->variance;
  const T varB = rhs.variance;

  Uncertain<T> result;

  result.mean = meanA + meanB;

  if (propagateUncertainty) {
    if (useCovariance) {
      covarianceMatrix.propagate(result.id,
				 this->id, 1,
				 rhs.id, 1);

      result.variance = covarianceMatrix.getCovariance(result.id, 
						       result.id);
    }
    else {
      result.variance = varA + varB;
    }
  }

  *this = result;

  return *this;
}

template <class T>
Uncertain<T> & 
Uncertain<T>::
operator+=(const T & rhs)
{
  *this += Uncertain<T>(rhs);

  return *this;
}

template <class T>
Uncertain<T> & 
Uncertain<T>::
operator-=(const Uncertain<T> & rhs)
{
  const T meanA = this->mean;
  const T meanB = rhs.mean;

  const T varA = this->variance;
  const T varB = rhs.variance;

  Uncertain<T> result;

  result.mean = meanA - meanB;

  if (propagateUncertainty) {
    if (useCovariance) {
      covarianceMatrix.propagate(result.id,
				 this->id, 1,
				 rhs.id, -1);

      result.variance = covarianceMatrix.getCovariance(result.id, 
						       result.id);
    }
    else {
      result.variance = varA + varB;
    }
  }

  *this = result;

  return *this;
}

template <class T>
Uncertain<T> & 
Uncertain<T>::
operator-=(const T & rhs)
{
  *this -= Uncertain<T>(rhs);

  return *this;
}

template <class T>
Uncertain<T> & 
Uncertain<T>::
operator*=(const Uncertain<T> & rhs)
{
  const T meanA = this->mean;
  const T meanB = rhs.mean;

  const T varA = this->variance;
  const T varB = rhs.variance;

  Uncertain<T> result;

  result.mean = meanA * meanB;

  if (propagateUncertainty) {
    T df_dA = meanB;
    T df_dB = meanA;

    if (useCovariance) {
      covarianceMatrix.propagate(result.id,
				 this->id, df_dA,
				 rhs.id, df_dB);

      result.variance = covarianceMatrix.getCovariance(result.id, 
						       result.id);
    }
    else {
      result.variance = 
	std::pow(meanA, 2) * varB +
	std::pow(meanB, 2) * varA +
	varA * varB;
    }
  }

  *this = result;

  return *this;
}

template <class T>
Uncertain<T> & 
Uncertain<T>::
operator*=(const T & rhs)
{
  *this *= Uncertain<T>(rhs);

  return *this;
}

template <class T>
Uncertain<T> & 
Uncertain<T>::
operator/=(const Uncertain<T> & rhs)
{
  const T meanA = this->mean;
  const T meanB = rhs.mean;

  const T varA = this->variance;
  const T varB = rhs.variance;

  Uncertain<T> result;

  result.mean = meanA / meanB;

  if (propagateUncertainty) {
    T df_dA = std::pow(meanB, -1);
    T df_dB = meanA * -1 * std::pow(meanB, -2);

    result.updateFirstOrder(this->id, varA, df_dA,
			    rhs.id, varB, df_dB);

    if (useSecondOrder) {
      T d2f_dA2  = 0;
      T d2f_dB2  = meanA * 2 * std::pow(meanB, -3);
      T d2f_dAdB = -1 * std::pow(meanB, -2);

      result.updateSecondOrder(varA, d2f_dA2,
			       d2f_dAdB,
			       varB, d2f_dB2);
    }
  }

  *this = result;

  return *this;
}

template <class T>
Uncertain<T> & 
Uncertain<T>::
operator/=(const T & rhs)
{
  *this /= Uncertain<T>(rhs);

  return *this;
}

/// Combines a Vector3 of means and a Vector3 of variances into a Vector3
/// of Uncertain values, each with a mean and variance.
///
template <class T>
Vector3<Uncertain<T> >
Uncertain<T>:: 
zip(Vector3<T> const & mean,
    Vector3<T> const & variance)
{
  Vector3<Uncertain<T> > zipped;

  for (int dim = 0; dim < 3; dim++) {
    zipped.set(dim, Uncertain<T>(mean.get(dim),
				 variance.get(dim)));
  }

  return zipped;
}

/// Combines a Matrix3x3 of means and a Matrix3x3 of variances into a Matrix3x3
/// of Uncertain values, each with a mean and variance.
///
template <class T>
Matrix3x3<Uncertain<T> > 
Uncertain<T>::
zip(Matrix3x3<T> const & mean,
    Matrix3x3<T> const & variance)
{
  Matrix3x3<Uncertain<T> > zipped;

  for (int component = 0; component < 3*3; component++) {
    zipped.set(component, Uncertain<T>(mean.get(component),
				       variance.get(component)));
  }

  return zipped;
}

/// Performs first-order propagation of variance for the function
/// x = f(A)
/// given the variance of A and the derivative df/dA.
///
template <class T>
T 
Uncertain<T>::
firstOrderVariance(T varA, T df_dA)
{
  T varX = varA * std::pow(df_dA, 2);

  return varX;
}

/// Performs first-order propagation of variance for the function
/// x = f(A, B)
/// given the variances of A and B and the partial derivatives df/dA and df/dB.
/// Assumes the covariance of A and B is 0.
///
template <class T>
T 
Uncertain<T>::
firstOrderVariance(T varA, T df_dA, 
		   T varB, T df_dB)
{
  T varX = 
    varA * std::pow(df_dA, 2) +
    varB * std::pow(df_dB, 2);

  return varX;
}

/// Performs first-order propagation of variance for the function
/// x = f(A, B)
/// given the variances of A and B, the partial derivatives df/dA and df/dB,
/// and the covariance of A and B.
///
template <class T>
T 
Uncertain<T>::
firstOrderVariance(T varA, T df_dA, 
		   T covAB,
		   T varB, T df_dB)
{
  T varX = 
    varA * std::pow(df_dA, 2) +
    2 * covAB * df_dA * df_dB +
    varB * std::pow(df_dB, 2);

  return varX;
}

/// Performs second-order propagation of mean for the function
/// x = f(A)
/// given the variance of A, the second derivative d^2f/dA^2,
/// and the first-order mean of x.
///
template <class T>
T Uncertain<T>::secondOrderMean(T firstOrderMean, 
				T varA, T d2f_dA2)
{
  T meanX = firstOrderMean + 0.5 * varA * d2f_dA2;

  return meanX;
}

/// Performs second-order propagation of mean for the function
/// x = f(A, B)
/// given the variances of A and B, the second partial derivatives 
/// d^2f/dA^2 and d^2f/dB^2,
/// and the first-order mean of x.
/// Assumes the covariance of A and B is 0.
///
template <class T>
T 
Uncertain<T>::
secondOrderMean(T firstOrderMean, 
		T varA, T d2f_dA2, 
		T varB, T d2f_dB2)
{
  T meanX = 
    firstOrderMean +
    0.5 * varA * d2f_dA2 +
    0.5 * varB * d2f_dB2;

  return meanX;
}

/// Performs second-order propagation of mean for the function
/// x = f(A, B)
/// given the variances of A and B, the second partial derivatives 
/// d^2f/dA^2, d^2f/dB^2 and d^2f/dAdB,
/// the covariance of A and B,
/// and the first-order mean of x.
///
template <class T>
T 
Uncertain<T>::
secondOrderMean(T firstOrderMean, 
		T varA, T d2f_dA2,
		T covAB, T d2f_dAdB, 
		T varB, T d2f_dB2)
{
  T meanX = 
    firstOrderMean +
    0.5 * varA * d2f_dA2 +
    covAB * d2f_dAdB +
    0.5 * varB * d2f_dB2;

  return meanX;
}

/// Performs second-order propagation of variance for the function
/// x = f(A)
/// given the variance of A, the second derivative d^2f/dA^2,
/// and the first-order variance of x.
///
template <class T>
T 
Uncertain<T>::
secondOrderVariance(T firstOrderVariance, 
		    T varA, T d2f_dA2)
{
  T varX = firstOrderVariance + 0.5 * std::pow(varA, 2) * std::pow(d2f_dA2, 2);

  return varX;
}

/// Performs second-order propagation of variance for the function
/// x = f(A, B)
/// given the variances of A and B, the second partial derivatives 
/// d^2f/dA^2, d^2f/dB^2, and d^2f/dAdB,
/// and the first-order variance of x.
/// Assumes the covariance of A and B is 0.
///
template <class T>
T 
Uncertain<T>::
secondOrderVariance(T firstOrderVariance, 
		    T varA, T d2f_dA2, 
		    T d2f_dAdB, 
		    T varB, T d2f_dB2)
{
  T varX =
    firstOrderVariance +
    0.5 * std::pow(varA, 2) * std::pow(d2f_dA2, 2) +
    (varA * varB) * std::pow(d2f_dAdB, 2) +
    0.5 * std::pow(varB, 2) * std::pow(d2f_dB2, 2);

  return varX;
}

/// Performs second-order propagation of variance for the function
/// x = f(A, B)
/// given the variances of A and B, the second partial derivatives 
/// d^2f/dA^2, d^2f/dB^2 and d^2f/dAdB,
/// the covariance of A and B,
/// and the first-order variance of x.
///
template <class T>
T 
Uncertain<T>::
secondOrderVariance(T firstOrderVariance, 
		    T varA, T d2f_dA2,
		    T covAB, T d2f_dAdB, 
		    T varB, T d2f_dB2)
{
  T varX =
    firstOrderVariance +
    0.5 * std::pow(varA, 2) * std::pow(d2f_dA2, 2) +
    2 * varA * covAB * d2f_dA2 * d2f_dAdB +
    std::pow(covAB, 2) * d2f_dA2 * d2f_dB2 +
    (varA * varB + std::pow(covAB, 2)) * std::pow(d2f_dAdB, 2) +
    2 * covAB * varB * d2f_dB2 * d2f_dAdB +
    0.5 * std::pow(varB, 2) * std::pow(d2f_dB2, 2);

  return varX;
}

template <class T>
void
Uncertain<T>::
updateFirstOrder(unsigned int idA, T varA, T df_dA)
{
  if (useCovariance) {
    covarianceMatrix.propagate(this->id,
			       idA, df_dA);

    this->variance = covarianceMatrix.getCovariance(this->id, 
						    this->id);
  }
  else {
    this->variance = firstOrderVariance(varA, df_dA);
  }
}

template <class T>
void
Uncertain<T>::
updateFirstOrder(unsigned int idA, T varA, T df_dA,
		 unsigned int idB, T varB, T df_dB)
{
  if (useCovariance) {
    covarianceMatrix.propagate(this->id,
			       idA, df_dA,
			       idB, df_dB);

    this->variance = covarianceMatrix.getCovariance(this->id, 
						    this->id);
  }
  else {
    this->variance = firstOrderVariance(varA, df_dA, 
					varB, df_dB);
  }
}

template <class T>
void
Uncertain<T>::
updateSecondOrder(T varA, T d2f_dA2)
{
  this->mean = secondOrderMean(this->mean, 
			       varA, d2f_dA2);
      
  this->variance = secondOrderVariance(this->variance, 
				       varA, d2f_dA2);
}

template <class T>
void
Uncertain<T>::
updateSecondOrder(T varA, T d2f_dA2,
		  T d2f_dAdB,
		  T varB, T d2f_dB2)
{
  this->mean = secondOrderMean(this->mean, 
			       varA, d2f_dA2,
			       varB, d2f_dB2);
      
  this->variance = secondOrderVariance(this->variance, 
				       varA, d2f_dA2,
				       d2f_dAdB, 
				       varB, d2f_dB2);
}

// ================================================================

template <class T>
Uncertain<T> fabs(const Uncertain<T> & A) 
{
  const T meanA = A.mean;
  const T varA  = A.variance;

  Uncertain<T> result;

  result.mean     = std::fabs(meanA);
  result.variance = 0;

  if (Uncertain<T>::propagateUncertainty) {
    T df_dA = (A.mean < 0) ? -1 : 1;

    result.updateFirstOrder(A.id, varA, df_dA);

    if (Uncertain<T>::useSecondOrder) {
      T d2f_dA2 = 0;

      result.updateSecondOrder(varA, d2f_dA2);
    }
  }

  return result;
}

template <class T>
Uncertain<T> log(const Uncertain<T> & A) 
{
  const T meanA = A.mean;
  const T varA  = A.variance;

  Uncertain<T> result;

  result.mean     = std::log(meanA);
  result.variance = 0;

  if (Uncertain<T>::propagateUncertainty) {
    T df_dA = 1 / meanA;

    result.updateFirstOrder(A.id, varA, df_dA);

    if (Uncertain<T>::useSecondOrder) {
      T d2f_dA2 = -std::pow(meanA, -2);

      result.updateSecondOrder(varA, d2f_dA2);
    }
  }

  return result;
}

template <class T>
Uncertain<T> exp(const Uncertain<T> & A) 
{
  const T meanA = A.mean;
  const T varA  = A.variance;

  Uncertain<T> result;

  result.mean     = std::exp(meanA);
  result.variance = 0;

  if (Uncertain<T>::propagateUncertainty) {
    T df_dA = std::exp(meanA);

    result.updateFirstOrder(A.id, varA, df_dA);

    if (Uncertain<T>::useSecondOrder) {
      T d2f_dA2 = std::exp(meanA);

      result.updateSecondOrder(varA, d2f_dA2);
    }
  }

  return result;
}

template <class T>
Uncertain<T> pow(const Uncertain<T> & A, const Uncertain<T> & B)
{
  const T meanA = A.mean;
  const T meanB = B.mean;

  const T varA = A.variance;
  const T varB = B.variance;

  //handle negative bases with integer exponents
  if (meanA < 0) {
    assert((int)meanB == meanB);

    if ((int)meanB % 2 == 0) {
      return pow(A * -1., B);
    }
    else {
      return pow(A * -1., B) * -1.;
    }
  }
  else {
    Uncertain<T> result;

    result.mean      = std::pow(meanA, meanB);
    result.variance  = 0;

    if (Uncertain<T>::propagateUncertainty) {
      T df_dA = meanB * std::pow(meanA, meanB - 1);
      T df_dB = std::log(meanA) * std::pow(meanA, meanB);

      if (meanA == 0) {
	df_dB = 0;
      }

      result.updateFirstOrder(A.id, varA, df_dA,
			      B.id, varB, df_dB);

      if (Uncertain<T>::useSecondOrder) {
	T d2f_dA2  = meanB * (meanB - 1) * std::pow(meanA, meanB - 2);
	T d2f_dB2  = std::pow(std::log(meanA), 2) * std::pow(meanA, meanB);
	T d2f_dAdB =
	  std::pow(meanA, meanB - 1) +
	  meanB * std::log(meanA) * std::pow(meanA, meanB - 1);

	if (meanA == 0) {
	  d2f_dB2  = 0;
	  d2f_dAdB = 0;
	}
	
	result.updateSecondOrder(varA, d2f_dA2,
				 d2f_dAdB,
				 varB, d2f_dB2);
      }
    }

    return result;
  }
}

template <class T>
Uncertain<T> pow(const Uncertain<T> & A, T b)
{
  return pow(A, Uncertain<T>(b));
}

template <class T>
Uncertain<T> pow(T a, const Uncertain<T> & B)
{
  return pow(Uncertain<T>(a), B);
}

template <class T>
Uncertain<T> sqrt(const Uncertain<T> & A) 
{
  const T meanA = A.mean;
  const T varA  = A.variance;

  Uncertain<T> result;

  result.mean     = std::sqrt(meanA);
  result.variance = 0;

  //can only take sqrt with mean == 0 if variance == 0, since otherwise 
  //sometimes the distribution would be negative
  if (meanA == 0) {
    assert(varA == 0);
  }
  else {
    if (Uncertain<T>::propagateUncertainty) {
      T df_dA = 0.5 * std::pow(meanA, -0.5);

      result.updateFirstOrder(A.id, varA, df_dA);

      if (Uncertain<T>::useSecondOrder) {
	T d2f_dA2 = -0.25 * std::pow(meanA, -1.5);

	result.updateSecondOrder(varA, d2f_dA2);
      }
    }
  }

  return result;
}

template <class T>
Uncertain<T> sin(const Uncertain<T> & A) 
{
  const T meanA = A.mean;
  const T varA  = A.variance;

  Uncertain<T> result;

  result.mean     = sin(meanA);
  result.variance = 0;

  if (Uncertain<T>::propagateUncertainty) {
    T df_dA = cos(meanA);

    result.updateFirstOrder(A.id, varA, df_dA);

    if (Uncertain<T>::useSecondOrder) {
      T d2f_dA2 = -sin(meanA);

      result.updateSecondOrder(varA, d2f_dA2);
    }
  }

  return result;
}

template <class T>
Uncertain<T> cos(const Uncertain<T> & A) 
{
  const T meanA = A.mean;
  const T varA  = A.variance;

  Uncertain<T> result;

  result.mean     = cos(meanA);
  result.variance = 0;

  if (Uncertain<T>::propagateUncertainty) {
    T df_dA = -sin(meanA);

    result.updateFirstOrder(A.id, varA, df_dA);

    if (Uncertain<T>::useSecondOrder) {
      T d2f_dA2 = -cos(meanA);

      result.updateSecondOrder(varA, d2f_dA2);
    }
  }

  return result;
}

template <class T>
Uncertain<T> atan2(const Uncertain<T> & A, const Uncertain<T> & B)
{
  const T meanA = A.mean;
  const T meanB = B.mean;

  const T varA = A.variance;
  const T varB = B.variance;

  Uncertain<T> result;

  result.mean      = atan2(meanA, meanB);
  result.variance  = 0;

  if (Uncertain<T>::propagateUncertainty) {
    const T B2A2 = std::pow(meanB, 2) + std::pow(meanA, 2);

    T df_dA = meanB / B2A2;
    T df_dB = -meanA / B2A2;

    result.updateFirstOrder(A.id, varA, df_dA,
			    B.id, varB, df_dB);

    if (Uncertain<T>::useSecondOrder) {
      T d2f_dA2  = meanB * -1 * std::pow(B2A2, -2) * 2 * meanA;
      T d2f_dB2  = -meanA * -1 * std::pow(B2A2, -2) * 2 * meanB;
      T d2f_dAdB =
	std::pow(B2A2, -1) + meanB * -1 * std::pow(B2A2, -2) * 2 * meanB;

      result.updateSecondOrder(varA, d2f_dA2,
			       d2f_dAdB,
			       varB, d2f_dB2);
    }
  }

  return result;
}

template <class T>
Uncertain<T> atan2(const Uncertain<T> & A, T b)
{
  return atan2(A, Uncertain<T>(b));
}

template <class T>
Uncertain<T> atan2(T a, const Uncertain<T> & B)
{
  return atan2(Uncertain<T>(a), B);
}

// ================================================================

template <class T>
std::ostream & operator<<(std::ostream & os, const Uncertain<T> & rhs)
{
  os << rhs.getMean() << " +/- " << rhs.getStdDev();

  return os;
}

template <class T>
bool operator==(const Uncertain<T> & A, const Uncertain<T> & B)
{
  return A.getMean() == B.getMean();
}

template <class T>
bool operator==(const Uncertain<T> & A, T b)
{
  return A == Uncertain<T>(b);
}

template <class T>
bool operator==(T a, const Uncertain<T> & B)
{
  return Uncertain<T>(a) == B;
}

template <class T>
bool operator!=(const Uncertain<T> & A, const Uncertain<T> & B)
{
  return !(A == B);
}

template <class T>
bool operator!=(const Uncertain<T> & A, T b)
{
  return A != Uncertain<T>(b);
}

template <class T>
bool operator!=(T a, const Uncertain<T> & B)
{
  return Uncertain<T>(a) != B;
}

template <class T>
bool operator<(const Uncertain<T> & A, const Uncertain<T> & B)
{
  return A.getMean() < B.getMean();
}

template <class T>
bool operator<(const Uncertain<T> & A, T b)
{
  return A < Uncertain<T>(b);
}

template <class T>
bool operator<(T a, const Uncertain<T> & B)
{
  return Uncertain<T>(a) < B;
}

template <class T>
bool operator>(const Uncertain<T> & A, const Uncertain<T> & B)
{
  return B < A;
}

template <class T>
bool operator>(const Uncertain<T> & A, T b)
{
  return A > Uncertain<T>(b);
}

template <class T>
bool operator>(T a, const Uncertain<T> & B)
{
  return Uncertain<T>(a) > B;
}

template <class T>
bool operator<=(const Uncertain<T> & A, const Uncertain<T> & B)
{
  return !(A > B);
}

template <class T>
bool operator<=(const Uncertain<T> & A, T b)
{
  return A <= Uncertain<T>(b);
}

template <class T>
bool operator<=(T a, const Uncertain<T> & B)
{
  return Uncertain<T>(a) <= B;
}

template <class T>
bool operator>=(const Uncertain<T> & A, const Uncertain<T> & B)
{
  return !(A < B);
}

template <class T>
bool operator>=(const Uncertain<T> & A, T b)
{
  return A >= Uncertain<T>(b);
}

template <class T>
bool operator>=(T a, const Uncertain<T> & B)
{
  return Uncertain<T>(a) >= B;
}

template <class T>
inline Uncertain<T> operator+(Uncertain<T> lhs, const Uncertain<T> & rhs)
{
  lhs += rhs;
  return lhs;
}

template <class T>
inline Uncertain<T> operator+(Uncertain<T> lhs, const T & rhs) 
{
  lhs += rhs;
  return lhs;
}

template <class T>
inline Uncertain<T> operator+(const T & lhs, const Uncertain<T> & rhs)
{
  Uncertain<T> uncertainLHS(lhs);

  uncertainLHS += rhs;

  return uncertainLHS;
}

template <class T>
inline Uncertain<T> operator-(Uncertain<T> lhs, const Uncertain<T> & rhs)
{
  lhs -= rhs;
  return lhs;
}

template <class T>
inline Uncertain<T> operator-(Uncertain<T> lhs, const T & rhs)
{
  lhs -= rhs;
  return lhs;
}

template <class T>
inline Uncertain<T> operator-(const T & lhs, const Uncertain<T> & rhs)
{
  Uncertain<T> uncertainLHS(lhs);

  uncertainLHS -= rhs;

  return uncertainLHS;
}

template <class T>
inline Uncertain<T> operator*(Uncertain<T> lhs, const Uncertain<T> & rhs)
{
  lhs *= rhs;
  return lhs;
}

template <class T>
inline Uncertain<T> operator*(Uncertain<T> lhs, const T & rhs)
{
  lhs *= rhs;
  return lhs;
}

template <class T>
inline Uncertain<T> operator*(const T & lhs, const Uncertain<T> & rhs)
{
  Uncertain<T> uncertainLHS(lhs);

  uncertainLHS *= rhs;

  return uncertainLHS;
}

template <class T>
inline Uncertain<T> operator/(Uncertain<T> lhs, const Uncertain<T> & rhs)
{
  lhs /= rhs;
  return lhs;
}

template <class T>
inline Uncertain<T> operator/(Uncertain<T> lhs, const T & rhs)
{
  lhs /= rhs;
  return lhs;
}

template <class T>
inline Uncertain<T> operator/(const T & lhs, const Uncertain<T> & rhs)
{
  Uncertain<T> uncertainLHS(lhs);

  uncertainLHS /= rhs;

  return uncertainLHS;
}

template <class T>
void Uncertain<T>::setCovariance(Uncertain<T> const & A, Uncertain<T> const & B, double covariance){
    Uncertain<T>::covarianceMatrix.setCovariance(A.id, B.id, covariance);
}

}

// ================================================================

#endif  // #ifndef UNCERTAIN_H_

