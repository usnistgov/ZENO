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
// Authors: Andrew Schultz <ajs42@buffalo.edu>
// Created: 2023-02-12
//
// ================================================================

#include "FfParser.h"

#include <boost/spirit/include/qi.hpp>

// ================================================================

using namespace zeno;

ff_parser::FfParser::FfParser
(std::istream &in,
 int activeModel,
 Potential<double> * potential) :
  in(in),
  potential(potential),
  activeModel(activeModel),
  ffSection(None) {

}

int ff_parser::FfParser::parse() {
  // Copy input file contents into string
  std::string storage;
  in.unsetf(std::ios::skipws);
  std::copy(
    std::istream_iterator<char>(in),
    std::istream_iterator<char>(),
    std::back_inserter(storage));
  
  std::string::const_iterator first = storage.begin();
  std::string::const_iterator last  = storage.end();

  bool parseResult = boost::spirit::qi::parse(first, last, FfParserGrammar(this));

  if (parseResult == false || first != last) {
    std::cerr << "Error parsing ff file: Unparseable: " << std::string(first, last) << std::endl;

    exit(1);
  }

  return 0;
}

void
ff_parser::FfParser::transformToLower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
        [](unsigned char c){ return std::tolower(c); });
}

void
ff_parser::FfParser::setSection(FfSection s) {
  ffSection = s;
}

void ff_parser::FfParser::processThreeInts(int a, int b, int c) {

  if (ffSection == Bonds) {
    potential->addBondPair(a, activeModel, b, c);
  }
  else {
    std::cout << "Found three ints outside bonds section" << std::endl;
  }
}

void ff_parser::FfParser::processFourInts(int a, int b, int c, int d) {

  if (ffSection == Angles) {
    potential->addAngleTriplet(a, activeModel, b, c, d);
  }
  else {
    std::cout << "Found four ints outside angles section" << std::endl;
  }
}


void ff_parser::FfParser::setBondStyleName(std::string bondStyleName) {
  transformToLower(bondStyleName);
  potential->setBondStyleName(bondStyleName);
}

void ff_parser::FfParser::setAngleStyleName(std::string angleStyleName) {
  transformToLower(angleStyleName);
  potential->setAngleStyleName(angleStyleName);
}

void ff_parser::FfParser::setBondCoeff4(int bondType, double c1, double c2, double c3, double c4) {
  double* coeffs = new double[4];
  coeffs[0] = c1;
  coeffs[1] = c2;
  coeffs[2] = c3;
  coeffs[3] = c4;
  potential->setBondCoeff(bondType, coeffs);
}

void ff_parser::FfParser::setBondCoeff2(int bondType, double c1, double c2) {
  double* coeffs = new double[2];
  coeffs[0] = c1;
  coeffs[1] = c2;
  potential->setBondCoeff(bondType, coeffs);
}

void ff_parser::FfParser::setAngleCoeff2(int angleType, double c1, double c2) {
  double* coeffs = new double[2];
  coeffs[0] = c1;
  coeffs[1] = c2;
  potential->setAngleCoeff(angleType, coeffs);
}

void ff_parser::FfParser::setNonbondStyleName(std::string nonbondStyleName) {
  transformToLower(nonbondStyleName);
  potential->setNonbondStyleName(nonbondStyleName);
}

void ff_parser::FfParser::setNonbondCoeff1(int sphereType, double c1) {
  double* coeffs = new double[1];
  coeffs[0] = c1;
  potential->setNonbondCoeff(sphereType, coeffs);
}

void ff_parser::FfParser::setNonbondCoeff2(int sphereType, double c1, double c2) {
  double* coeffs = new double[2];
  coeffs[0] = c1;
  coeffs[1] = c2;
  potential->setNonbondCoeff(sphereType, coeffs);
}
