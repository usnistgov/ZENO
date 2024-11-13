// ================================================================
//
// This software is not subject to copyright protection in the United States
// and is considered to be in the public domain. Permission to freely use,
// copy, modify, and distribute this software and its documentation without fee
// is hereby granted, provided that this notice and disclaimer of warranty
// appears in all copies.
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
// Authors: Andrew Schultz <ajs42@buffalo.edu>
// Created: Wed Feb 22 2023
//
// ================================================================

#include "Potential.h"

using namespace zeno;

template <class T>
Potential<T>::Potential()
: bondStyle(BondStyle::None),
  nonbondStyle(NonbondStyle::HardSphere) {

}

template <class T>
Potential<T>::Potential(const Potential &p) : bondStyle(p.getBondStyle()),
                                           nonbondStyle(p.getNonbondStyle()),
                                           bondPairs(p.getBondPairs()) {
  for (int i=0; i<(int)p.getBondCoeffs()->size(); i++) {
    double* c = p.getBondCoeffs()->at(i);
    if (!c) {
      bondCoeffs.push_back(c);
    }
    else {
      if (bondStyle == BondStyle::Harmonic) {
        double* cNew = new double[2];
        cNew[0] = c[0];
        cNew[1] = c[1];
        bondCoeffs.push_back(cNew);
      }
    }
  }

  for (int i=0; i<(int)p.getNonbondCoeffs()->size(); i++) {
    double* c = p.getNonbondCoeffs()->at(i);
    if (!c) {
      nonbondCoeffs1.push_back(c);
    }
    else {
      if (nonbondStyle == NonbondStyle::HardSphere) {
        double* cNew = new double[1];
        cNew[0] = c[0];
        nonbondCoeffs1.push_back(cNew);
      }
      if (nonbondStyle == NonbondStyle::LennardJones) {
        double* cNew = new double[2];
        cNew[0] = c[0];
        cNew[1] = c[1];
        nonbondCoeffs1.push_back(cNew);
      }
    }
  }
}

template <class T>
Potential<T>::~Potential() {
  for (int i=0; i<(int)bondCoeffs.size(); i++) delete [] bondCoeffs[i];
}

/// Add a bond to the model.
///
template <class T>
void
Potential<T>::addBondPair(int bondType, int sphere1, int sphere2) {
  struct BondedPair bp = {bondType, sphere1, sphere2};
  bondPairs.push_back(bp);
}

/// Set the bond potential parameters for the given bond type.
///
template <class T>
void
Potential<T>::setBondCoeff2(int bondType, double c1, double c2) {
  if ((int)bondCoeffs.size() <= bondType) {
    bondCoeffs.resize(bondType+1, nullptr);
  }
  double* coeffs = new double[2];
  coeffs[0] = c1;
  coeffs[1] = c2;
  bondCoeffs[bondType] = coeffs;
}

/// Set the bond potential parameters for the given bond type.
///
template <class T>
void
Potential<T>::setNonbondCoeff2(int sphereType, double c1, double c2) {
  if ((int)nonbondCoeffs.size() <= sphereType) {
    nonbondCoeffs1.resize(sphereType+1, nullptr);
  }
  double* coeffs = new double[2];
  coeffs[0] = c1;
  coeffs[1] = c2;
  nonbondCoeffs1[sphereType] = coeffs;
}

/// Set the bond style for the model.
///
template <class T>
void
Potential<T>::setBondStyleName(std::string bondStyleName) {
  if (bondStyleName.compare("harmonic") == 0) {
    bondStyle = BondStyle::Harmonic;
  }
  else if (bondStyleName.compare("fene") == 0) {
    bondStyle = BondStyle::FENE;
  }
  else {
    std::cerr << "Unrecognized bond style " << bondStyleName << std::endl;

    exit(1);
  }
}

/// Set the bond style for the model.
///
template <class T>
void
Potential<T>::setNonbondStyleName(std::string nonbondStyleName) {
  if (nonbondStyleName.compare("hs") == 0) {
    nonbondStyle = NonbondStyle::HardSphere;
  }
  else if (nonbondStyleName.compare("lj") == 0) {
    nonbondStyle = NonbondStyle::LennardJones;
  }
  else {
    std::cerr << "Unrecognized nonbond style " << nonbondStyleName << std::endl;

    exit(1);
  }
}

template <class T>
const BondStyle
Potential<T>::getBondStyle() const {
  return bondStyle;
}

template <class T>
const NonbondStyle
Potential<T>::getNonbondStyle() const {
  return nonbondStyle;
}

template <class T>
std::vector<BondedPair>
Potential<T>::getBondPairs() const {
  return bondPairs;
}

template <class T>
const std::vector<double*> *
Potential<T>::getBondCoeffs() const {
  return &bondCoeffs;
}

template <class T>
const std::vector<double*> *
Potential<T>::getNonbondCoeffs() const {
  return &nonbondCoeffs1;
}

template <class T>
void
Potential<T>::initialize() {

  std::vector<PairPotential*> bondPotentials;

  for (BondedPair bp : bondPairs) {
    if ((int)bondPotentials.size() <= bp.bondType) {
      bondPotentials.resize(bp.bondType+1);
    }
    PairPotential* p2 = nullptr;
    if (bondPotentials[bp.bondType] == nullptr) {
      if (bondStyle == BondStyle::Harmonic) {
        double* iBondCoeff = bondCoeffs[bp.bondType];
        bondPotentials[bp.bondType] = new PairPotentialHarmonic(iBondCoeff[0], iBondCoeff[1]);
      }
    }
    else {
      p2 = bondPotentials[bp.bondType];
    }
    BondedPartner partner1 = {bp.sphere2, p2};
    bondedPartners[bp.sphere1].push_back(partner1);
    BondedPartner partner2 = {bp.sphere1, p2};
    bondedPartners[bp.sphere2].push_back(partner2);
  }

}

PairPotentialHarmonic::PairPotentialHarmonic(double r0_, double k_) : r0(r0_), k(k_) {
}

double PairPotentialHarmonic::u(double r2) {
  double dr = std::sqrt(r2) - r0;
  return k*dr*dr;
}

