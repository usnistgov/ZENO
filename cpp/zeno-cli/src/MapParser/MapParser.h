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

#ifndef map_parserMapParser_h_included
#define map_parserMapParser_h_included

#include <string>
#include <unordered_map>
#include <vector>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/qi_omit.hpp>
#include <boost/spirit/include/phoenix.hpp>

namespace map_parser {
    namespace qi = boost::spirit::qi;
    namespace phx = boost::phoenix;

    // Begin class MapParser
    class MapParser {
        public:
            MapParser(std::istream & in,
		              std::unordered_map<std::string, double> * atomIdToRadius);
          
            int parse();

        private:
            void addMapping(std::string atomId, double radius);
     
            std::istream & in;

            std::unordered_map<std::string, double> * atomIdToRadius;

            // Begin struct Mapping
            struct Mapping {
                std::string atomId;
                double radius;

                Mapping() {}
                Mapping(std::vector<char> atomIdVec, double radius)
                    : atomId(atomIdVec.begin(), atomIdVec.end()), radius(radius) {}
            };
            // End struct Mapping

            // Begin class MapParserGrammar
            class MapParserGrammar
                : public qi::grammar<std::string::const_iterator, std::vector<Mapping>()> {
    
                public:
                    MapParserGrammar() : MapParserGrammar::base_type(start) {
                        comment = (qi::lit('#') >> *(qi::char_ - qi::eol) >> qi::omit[+qi::eol])
                                  [qi::_val = qi::_1];

                        mapping = (+(qi::char_ - qi::blank) >> qi::omit[+qi::blank] >> qi::double_ >> qi::omit[*qi::space])
                                  [qi::_val = phx::construct<Mapping>(qi::_1, qi::_2)];

                        start = *(
                                  comment
                                  |
                                  mapping [phx::push_back(qi::_val, qi::_1)]
                                 );
                    }

                    qi::rule<std::string::const_iterator, std::vector<char>()> comment;
                    qi::rule<std::string::const_iterator, Mapping()> mapping;
                    qi::rule<std::string::const_iterator, std::vector<Mapping>()> start;
            };
            // End class MapParserGrammar
    };
    // End class MapParser
}

#endif
