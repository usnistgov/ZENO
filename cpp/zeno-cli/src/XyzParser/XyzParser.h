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

#ifndef xyz_parserXyzParser_h_included
#define xyz_parserXyzParser_h_included

#include <string>
#include <unordered_map>
#include <list>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/qi_omit.hpp>
#include <boost/spirit/include/phoenix.hpp>

#include "Geometry/MixedModel.h"

namespace xyz_parser {
    namespace qi = boost::spirit::qi;
    namespace phx = boost::phoenix;

    // Begin class XyzParser
    class XyzParser {            
        public:
            XyzParser(std::unordered_map<std::string, double> const & atomIdToRadius,
                      std::istream & in,
                      std::list<zeno::MixedModel<double> > * snapshots,
	                  std::list<std::string> * comments);
	    
            int parse();

        private:
	        void addSnapshot(int numAtoms);
	        void addComment(std::string comment);
	        void addAtom(std::string id,
		                 double x, double y, double z);

            std::istream & in;

            std::unordered_map<std::string, double> const * atomIdToRadius;
            std::list<zeno::MixedModel<double> > * snapshots;
            std::list<std::string> * comments;
            
            int numAtomsExpected;
            int numAtomsSeen;

            // Begin struct Atom
            struct Atom {
                std::string id;
                double x, y, z;

                Atom() {}
                Atom(std::vector<char> idVec, double x, double y, double z)
                    : id(idVec.begin(), idVec.end()), x(x), y(y), z(z) {}
            };
            // End struct Atom

            // Begin struct Snapshot
            struct Snapshot {
                int numAtoms;
                std::string comment;
                std::vector<Atom> atoms;

                Snapshot() {}
                Snapshot(int numAtoms, std::vector<char> commentVec, std::vector<Atom> atoms)
                    : numAtoms(numAtoms), comment(commentVec.begin(), commentVec.end()), atoms(atoms) {}
            };
            // End struct Snapshot

            // Begin class XyzParserGrammar
            class XyzParserGrammar 
                : public qi::grammar<std::string::const_iterator, std::vector<Snapshot>()> {

                public:
                    XyzParserGrammar() : XyzParserGrammar::base_type(start) {
                        comment = (*(qi::char_ - qi::eol) >> +qi::eol)
                                  [qi::_val = qi::_1];

                        atom = (+(qi::char_ - qi::blank) >> qi::omit[+qi::blank] >>
                                qi::double_ >> qi::omit[+qi::blank] >>
                                qi::double_ >> qi::omit[+qi::blank] >>
                                qi::double_ >> qi::omit[*qi::space])
                               [qi::_val = phx::construct<Atom>(qi::_1, qi::_2, qi::_3, qi::_4)];

                        snapshot = (qi::int_ >> qi::omit[*qi::space] >> 
                                    comment >>
                                    *(atom))
                                   [qi::_val = phx::construct<Snapshot>(qi::_1, qi::_2, qi::_3)];

                        start = (*(snapshot))
                                [qi::_val = qi::_1];
                    }

                    qi::rule<std::string::const_iterator, std::vector<char>()> comment;
                    qi::rule<std::string::const_iterator, Atom()> atom;
                    qi::rule<std::string::const_iterator, Snapshot()> snapshot;
                    qi::rule<std::string::const_iterator, std::vector<Snapshot>()> start;
            };
            // End class XyzParserGrammar
    };
    // End class XyzParser
}

#endif
