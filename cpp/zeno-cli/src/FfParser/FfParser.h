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
// Created: 2023-02-08
//
// ================================================================

#ifndef ff_parserFfParser_h_included
#define ff_parserFfParser_h_included

#include <vector>
#include <string>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/spirit/include/qi_no_case.hpp>
#include <boost/spirit/include/qi_omit.hpp>
#include <boost/phoenix/bind/bind_member_function.hpp>

#include "Potential.h"

namespace ff_parser {
    namespace qi = boost::spirit::qi;
    namespace phx = boost::phoenix;

    // Begin class FfParser
    class FfParser {
        public:
            enum FfSection {None, Bonds, Angles, Dihedrals};
            FfParser(std::istream &in, int activeModel, zeno::Potential<double> * potential);

            int parse();

        private:
            void transformToLower(std::string s);
	    void setSection(FfSection s);
	    void processThreeInts(int a, int b, int c);
	    void processFourInts(int a, int b, int c, int d);
	    void setBondStyleName(std::string bondStyleName);
	    void setAngleStyleName(std::string angleStyleName);
	    void setBondCoeff2(int bondType, double c1, double c2);
	    void setAngleCoeff2(int bondType, double c1, double c2);
	    void setBondCoeff4(int bondType, double c1, double c2, double c3, double c4);
	    void setNonbondStyleName(std::string nonbondStyleName);
	    void setNonbondCoeff1(int sphereType, double c1);
	    void setNonbondCoeff2(int sphereType, double c1, double c2);

            std::istream & in;
            zeno::Potential<double> * potential;
            int activeModel;
            FfSection ffSection;

            // Begin class FfParserGrammar
            class FfParserGrammar
                : public qi::grammar<std::string::const_iterator, void()> {

                public:
                    FfParserGrammar(ff_parser::FfParser *parent) : FfParserGrammar::base_type(start) {

                        word %= +(qi::char_ - qi::space);

                        bondsSection = (qi::no_case["bonds"] >> qi::omit[*qi::blank] >> qi::eol)
                                 [phx::bind(&ff_parser::FfParser::setSection, parent,
                                     FfSection::Bonds)];

                        anglesSection = (qi::no_case["angles"] >> qi::omit[*qi::blank] >> qi::eol)
                                 [phx::bind(&ff_parser::FfParser::setSection, parent,
                                     FfSection::Angles)];

                        dihedralsSection = (qi::no_case["dihedrals"] >> qi::omit[*qi::blank] >> qi::eol)
                                 [phx::bind(&ff_parser::FfParser::setSection, parent,
                                     FfSection::Dihedrals)];

                        fourInts = (qi::int_ >> qi::omit[+qi::blank] >>
                                    qi::int_ >> qi::omit[+qi::blank] >>
                                    qi::int_ >> qi::omit[+qi::blank] >>
                                    qi::int_ >> qi::eol)
                                 [phx::bind(&ff_parser::FfParser::processFourInts, parent, 
                                            qi::_1, qi::_2, qi::_3, qi::_4)];

                        threeInts = (qi::int_ >> qi::omit[+qi::blank] >>
                                     qi::int_ >> qi::omit[+qi::blank] >>
                                     qi::int_ >> qi::eol)
                                 [phx::bind(&ff_parser::FfParser::processThreeInts, parent, 
                                            qi::_1, qi::_2, qi::_3)];

                        bondStyle = (qi::no_case["bond_style"] >> qi::omit[+qi::blank] >>
                                  word >> qi::omit[*qi::space])
                                 [phx::bind(&ff_parser::FfParser::setBondStyleName, parent, 
                                            qi::_1)];

                        angleStyle = (qi::no_case["angle_style"] >> qi::omit[+qi::blank] >>
                                  word >> qi::omit[*qi::space])
                                 [phx::bind(&ff_parser::FfParser::setAngleStyleName, parent, 
                                            qi::_1)];

                        bondCoeff2 = (qi::no_case["bond_coeff"] >> qi::omit[+qi::blank] >>
                                  qi::int_ >> qi::omit[+qi::blank] >>
                                  qi::double_ >> qi::omit[+qi::blank] >>
                                  qi::double_ >> qi::omit[*qi::space])
                                 [phx::bind(&ff_parser::FfParser::setBondCoeff2, parent, 
                                            qi::_1, qi::_2, qi::_3)];

                        bondCoeff4 = (qi::no_case["bond_coeff"] >> qi::omit[+qi::blank] >>
                                  qi::int_ >> qi::omit[+qi::blank] >>
                                  qi::double_ >> qi::omit[+qi::blank] >>
                                  qi::double_ >> qi::omit[+qi::blank] >>
                                  qi::double_ >> qi::omit[+qi::blank] >>
                                  qi::double_ >> qi::omit[*qi::space])
                                 [phx::bind(&ff_parser::FfParser::setBondCoeff4, parent, 
                                            qi::_1, qi::_2, qi::_3, qi::_4, qi::_5)];

                        angleCoeff2 = (qi::no_case["angle_coeff"] >> qi::omit[+qi::blank] >>
                                  qi::int_ >> qi::omit[+qi::blank] >>
                                  qi::double_ >> qi::omit[+qi::blank] >>
                                  qi::double_ >> qi::omit[*qi::space])
                                 [phx::bind(&ff_parser::FfParser::setAngleCoeff2, parent, 
                                            qi::_1, qi::_2, qi::_3)];

                        nonbondStyle = (qi::no_case["nonbond_style"] >> qi::omit[+qi::blank] >>
                                  word >> qi::omit[*qi::space])
                                 [phx::bind(&ff_parser::FfParser::setNonbondStyleName, parent, 
                                            qi::_1)];

                        nonbondCoeff2 = (qi::no_case["nonbond_coeff"] >> qi::omit[+qi::blank] >>
                                  qi::int_ >> qi::omit[+qi::blank] >>
                                  qi::double_ >> qi::omit[+qi::blank] >>
                                  qi::double_ >> qi::omit[*qi::space])
                                 [phx::bind(&ff_parser::FfParser::setNonbondCoeff2, parent, 
                                            qi::_1, qi::_2, qi::_3)];

                        nonbondCoeff1 = (qi::no_case["nonbond_coeff"] >> qi::omit[+qi::blank] >>
                                  qi::int_ >> qi::omit[+qi::blank] >>
                                  qi::double_ >> qi::omit[*qi::space])
                                 [phx::bind(&ff_parser::FfParser::setNonbondCoeff1, parent, 
                                            qi::_1, qi::_2)];

                        section = bondsSection |
                                  anglesSection |
                                  dihedralsSection;

                        sectionContent = fourInts |
                                         threeInts;

                        parameter = bondStyle |
                                    angleStyle |
                                    bondCoeff4 |
                                    bondCoeff2 |
                                    angleCoeff2 |
                                    nonbondStyle |
                                    nonbondCoeff2 |
                                    nonbondCoeff1;

                        start = *(section |
                                  sectionContent |
                                  parameter);

                    }

                    qi::rule<std::string::const_iterator, std::string()> word;
                    qi::rule<std::string::const_iterator, void()> section;
                    qi::rule<std::string::const_iterator, void()> bondsSection;
                    qi::rule<std::string::const_iterator, void()> anglesSection;
                    qi::rule<std::string::const_iterator, void()> dihedralsSection;
                    qi::rule<std::string::const_iterator, void()> sectionContent;
                    qi::rule<std::string::const_iterator, void()> threeInts;
                    qi::rule<std::string::const_iterator, void()> fourInts;
                    qi::rule<std::string::const_iterator, void()> bondStyle;
                    qi::rule<std::string::const_iterator, void()> angleStyle;
                    qi::rule<std::string::const_iterator, void()> bondCoeff4;
                    qi::rule<std::string::const_iterator, void()> bondCoeff2;
                    qi::rule<std::string::const_iterator, void()> angleCoeff2;
                    qi::rule<std::string::const_iterator, void()> nonbondStyle;
                    qi::rule<std::string::const_iterator, void()> nonbondCoeff2;
                    qi::rule<std::string::const_iterator, void()> nonbondCoeff1;
                    qi::rule<std::string::const_iterator, void()> parameter;
                    qi::rule<std::string::const_iterator, void()> start;
            };
            // End class FfParserGrammar
    };
    // End class FfParser
}

#endif
