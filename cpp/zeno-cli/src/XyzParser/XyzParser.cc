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
// Created: 2023-01-18
//
// ================================================================

#include "XyzParser.h"

#include <boost/spirit/include/qi.hpp>

#include <iostream>

// ================================================================

using namespace zeno;

xyz_parser::XyzParser::XyzParser
(std::unordered_map<std::string, double> const & atomIdToRadius,
 std::istream & in,
 std::list<zeno::MixedModel<double> > * snapshots,
 std::list<std::string> * comments)
  : in(in),
    atomIdToRadius(&atomIdToRadius),
    snapshots(snapshots),
    comments(comments),
    numAtomsExpected(0),
    numAtomsSeen(0) {

}

int xyz_parser::XyzParser::parse() {
  // Copy input file contents into string
  std::string storage;
  in.unsetf(std::ios::skipws);
  std::copy(
    std::istream_iterator<char>(in),
    std::istream_iterator<char>(),
    std::back_inserter(storage));
  
  std::string::const_iterator first = storage.begin();
  std::string::const_iterator last  = storage.end();

  std::vector<Snapshot> parsedSnapshots;

  bool parseResult = boost::spirit::qi::parse(first, last, XyzParserGrammar(), parsedSnapshots);

  if (parseResult == false || first != last) {
    std::cerr << "Error parsing XYZ file: Unparseable: " << std::string(first, last) << std::endl;

    exit(1);
  }

  for (auto snapshot = parsedSnapshots.begin(); snapshot != parsedSnapshots.end(); ++snapshot) {
    this->addSnapshot(snapshot->numAtoms);

    this->addComment(snapshot->comment);

    for (auto atom = snapshot->atoms.begin(); atom != snapshot->atoms.end(); ++atom) {
      this->addAtom(atom->id, atom->x, atom->y, atom->z);
    }
  }

  return 0;
}

void xyz_parser::XyzParser::addSnapshot(int numAtoms) {
  if (numAtomsSeen < numAtomsExpected) {
    std::cerr << "Error parsing XYZ file: Expected " << numAtomsExpected
	      << " atoms in snapshot but found only " << numAtomsSeen
	      << std::endl;

    exit(1);
  }

  numAtomsExpected = numAtoms;

  numAtomsSeen = 0;

  snapshots->emplace_back();

  comments->emplace_back();
}

void xyz_parser::XyzParser::addComment(std::string comment) {
  comments->back() = comment;
}

void xyz_parser::XyzParser::addAtom(std::string id,
				                    double x, double y, double z) {

  ++ numAtomsSeen;

  if (numAtomsSeen > numAtomsExpected) {
    std::cerr << "Error parsing XYZ file: Expected only " << numAtomsExpected
	      << " atoms in snapshot but found more"
	      << std::endl;

    exit(1);
  }

  double radius{};
  
  try {
    radius = atomIdToRadius->at(id);
  }
  catch (std::out_of_range const & oor) {
    std::cerr << "Error parsing XYZ file: Found unknown atom ID " << id
	      << std::endl;

    exit(1);
  }

  Vector3<double> center(x, y, z);

  Sphere<double> sphere(center, radius);

  snapshots->back().addSphere(sphere);
}
