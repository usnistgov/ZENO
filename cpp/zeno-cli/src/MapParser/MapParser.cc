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
// Created: 2023-01-13
//
// ================================================================

#include "MapParser.h"

#include <boost/spirit/include/qi.hpp>

#include <iostream>

// ================================================================

map_parser::MapParser::MapParser
(std::istream & in,
 std::unordered_map<std::string, double> * atomIdToRadius)
  : in(in),
    atomIdToRadius(atomIdToRadius) {

}

int map_parser::MapParser::parse() {  
  // Copy input file contents into string
  std::string storage;
  in.unsetf(std::ios::skipws);
  std::copy(
    std::istream_iterator<char>(in),
    std::istream_iterator<char>(),
    std::back_inserter(storage));
  
  std::string::const_iterator first = storage.begin();
  std::string::const_iterator last  = storage.end();

  std::vector<Mapping> mappings;

  bool parseResult = boost::spirit::qi::parse(first, last, MapParserGrammar(), mappings);

  if (parseResult == false || first != last) {
    std::cerr << "Error parsing MAP file: Unparseable: " << std::string(first, last) << std::endl;

    exit(1);
  }

  for (auto mapping = mappings.begin(); mapping != mappings.end(); ++mapping) {
    addMapping(mapping->atomId, mapping->radius);
  }

  return 0;
}

void map_parser::MapParser::addMapping(std::string atomId, double radius) {
  // Ignore atoms of radius 0
  if (radius != 0) {
    auto emplaceResult = atomIdToRadius->emplace(atomId, radius);

    if (emplaceResult.second == false) {
      std::cerr << "Error parsing MAP file: Found duplicate atom ID " << atomId
		<< " with radius " << radius << std::endl;

      exit(1);
    }
  }
}
