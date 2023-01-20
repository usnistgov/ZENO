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
// Created: 2023-01-19
//
// ================================================================

#ifndef bod_parserBodParser_h_included
#define bod_parserBodParser_h_included

#include <vector>
#include <string>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/spirit/include/qi_no_case.hpp>
#include <boost/spirit/include/qi_omit.hpp>
#include <boost/phoenix/bind/bind_member_function.hpp>

#include "ParametersWalkOnSpheres.h"
#include "ParametersInteriorSampling.h"
#include "ParametersResults.h"
#include "../ParametersLocal.h"
#include "Geometry/MixedModel.h"

namespace bod_parser {
    namespace qi = boost::spirit::qi;
    namespace phx = boost::phoenix;

    // Begin class BodParser
    class BodParser {
        public:
            BodParser(ParametersLocal * parametersLocal,
	                  std::istream &in,
	                  zeno::ParametersWalkOnSpheres * parametersWalkOnSpheres,
	                  zeno::ParametersInteriorSampling * parametersInteriorSampling,
	                  zeno::ParametersResults * parametersResults,
	                  zeno::MixedModel<double> * model);

            int parse();

        private:
            void addSphere(zeno::Vector3<double> center, double r);
	        void addCube(zeno::Vector3<double> minCoords, double s);
	        void addCuboid(zeno::Vector3<double> corner1,
		                   zeno::Vector3<double> corner2);
	        void addTriangle(zeno::Vector3<double> v1,
                             zeno::Vector3<double> v2,
                             zeno::Vector3<double> v3);
	        void addVoxels(std::string voxelsFileName);
	        void addTrajectory(std::string xyzFileName, std::string mapFileName);
	        void setST(double skinThickness);
	        void setRLAUNCH(double launchRadius);
	        void setHUNITS(double number, std::string unitString);
	        void setUNITS(std::string unitString);
	        void setTEMP(double number, std::string unitString);
	        void setMASS(double number, std::string unitString);
	        void setVISCOSITY(double number, std::string unitString);
	        void setBF(double buoyancyFactor);

            ParametersLocal * parametersLocal;
            std::istream & in;
            zeno::ParametersWalkOnSpheres * parametersWalkOnSpheres;
            zeno::ParametersInteriorSampling * parametersInteriorSampling;
            zeno::ParametersResults * parametersResults;
            zeno::MixedModel<double> * model;

            // Begin class BodParserGrammar
            class BodParserGrammar
                : public qi::grammar<std::string::const_iterator, void()> {

                public:
                    BodParserGrammar(bod_parser::BodParser *parent) : BodParserGrammar::base_type(start) {
                        vec3 = (qi::double_ >> qi::omit[+qi::blank] >> 
                                qi::double_ >> qi::omit[+qi::blank] >>
                                qi::double_)
                               [qi::_val = phx::construct<zeno::Vector3<double> >(qi::_1, qi::_2, qi::_3)];

                        word %= +(qi::char_ - qi::space);

                        sphere = (qi::no_case["SPHERE"] >> qi::omit[+qi::blank] >>
                                  vec3 >> qi::omit[+qi::blank] >>
                                  qi::double_ >> qi::omit[*qi::space])
                                 [phx::bind(&bod_parser::BodParser::addSphere, parent, 
                                            qi::_1, qi::_2)];

                        cube = (qi::no_case["CUBE"] >> qi::omit[+qi::blank] >>
                                vec3 >> qi::omit[+qi::blank] >> 
                                qi::double_ >> qi::omit[*qi::space])
                               [phx::bind(&bod_parser::BodParser::addCube, parent, 
                                          qi::_1, qi::_2)];

                        cuboid = (qi::no_case["CUBOID"] >> qi::omit[+qi::blank] >>
                                  vec3 >> qi::omit[+qi::blank] >> 
                                  vec3 >> qi::omit[*qi::space])
                                 [phx::bind(&bod_parser::BodParser::addCuboid, parent, 
                                            qi::_1, qi::_2)];

                        atom = (qi::no_case["ATOM"] >> qi::omit[+qi::blank] >>
                                qi::omit[qi::int_] >> qi::omit[+qi::blank] >>
                                qi::omit[word] >> qi::omit[+qi::blank] >>
                                qi::omit[word] >> qi::omit[+qi::blank] >>
                                -(qi::omit[word] >> qi::omit[+qi::blank]) >>
                                qi::omit[qi::int_] >> qi::omit[+qi::blank] >>
                                vec3 >> qi::omit[+qi::blank] >>
                                qi::omit[qi::double_] >> qi::omit[+qi::blank] >>
                                qi::double_ >> qi::omit[*qi::space])
                               [phx::bind(&bod_parser::BodParser::addSphere, parent, 
                                          qi::_1, qi::_2)];

                        voxels = (qi::no_case["VOXELS"] >> qi::omit[+qi::blank] >>
                                  word >> qi::omit[*qi::space])
                                 [phx::bind(&bod_parser::BodParser::addVoxels, parent,
                                  qi::_1)];

                        trajectory = (qi::no_case["TRAJECTORY"] >> qi::omit[+qi::blank] >>
                                      word >> qi::omit[+qi::blank] >> 
                                      word >> qi::omit[*qi::space])
                                     [phx::bind(&bod_parser::BodParser::addTrajectory, parent,
                                      qi::_1, qi::_2)];

                        parameter_st = (qi::no_case["ST"] >> qi::omit[+qi::blank] >>
                                        qi::double_ >> qi::omit[*qi::space])
                                       [phx::bind(&bod_parser::BodParser::setST, parent,
                                        qi::_1)];

                        parameter_rlaunch = (qi::no_case["RLAUNCH"] >> qi::omit[+qi::blank] >>
                                             qi::double_ >> qi::omit[*qi::space])
                                            [phx::bind(&bod_parser::BodParser::setRLAUNCH, parent,
                                             qi::_1)];

                        parameter_hunits = (qi::no_case["HUNITS"] >> qi::omit[+qi::blank] >>
                                            qi::double_ >> qi::omit[+qi::blank] >>
                                            word >> qi::omit[*qi::space])
                                           [phx::bind(&bod_parser::BodParser::setHUNITS, parent,
                                            qi::_1, qi::_2)];

                        parameter_units = (qi::no_case["UNITS"] >> qi::omit[+qi::blank] >>
                                           word >> qi::omit[*qi::space])
                                          [phx::bind(&bod_parser::BodParser::setUNITS, parent,
                                           qi::_1)];

                        parameter_temp = (qi::no_case["TEMP"] >> qi::omit[+qi::blank] >>
                                          qi::double_ >> qi::omit[+qi::blank] >>
                                          word >> qi::omit[*qi::space])
                                         [phx::bind(&bod_parser::BodParser::setTEMP, parent,
                                          qi::_1, qi::_2)];

                        parameter_mass = (qi::no_case["MASS"] >> qi::omit[+qi::blank] >>
                                          qi::double_ >> qi::omit[+qi::blank] >>
                                          word >> qi::omit[*qi::space])
                                         [phx::bind(&bod_parser::BodParser::setMASS, parent,
                                          qi::_1, qi::_2)];

                        parameter_viscosity = (qi::no_case["VISCOSITY"] >> qi::omit[+qi::blank] >>
                                               qi::double_ >> qi::omit[+qi::blank] >>
                                               word >> qi::omit[*qi::space])
                                              [phx::bind(&bod_parser::BodParser::setVISCOSITY, parent,
                                               qi::_1, qi::_2)];

                        parameter_bf = (qi::no_case["BF"] >> qi::omit[+qi::blank] >>
                                        qi::double_ >> qi::omit[*qi::space])
                                       [phx::bind(&bod_parser::BodParser::setBF, parent,
                                        qi::_1)];

                        parameter = parameter_st |
                                    parameter_rlaunch |
                                    parameter_hunits |
                                    parameter_units |
                                    parameter_temp |
                                    parameter_mass |
                                    parameter_viscosity |
                                    parameter_bf;

                        start = *(sphere |
                                  cube |
                                  cuboid |
                                  atom |
                                  voxels |
                                  trajectory |
                                  parameter);
                    }

                    qi::rule<std::string::const_iterator, zeno::Vector3<double>()> vec3;
                    qi::rule<std::string::const_iterator, std::string()> word;
                    qi::rule<std::string::const_iterator, void()> sphere;
                    qi::rule<std::string::const_iterator, void()> cube;
                    qi::rule<std::string::const_iterator, void()> cuboid;
                    qi::rule<std::string::const_iterator, void()> atom;
                    qi::rule<std::string::const_iterator, void()> voxels;
                    qi::rule<std::string::const_iterator, void()> trajectory;
                    qi::rule<std::string::const_iterator, void()> parameter_st;
                    qi::rule<std::string::const_iterator, void()> parameter_rlaunch;
                    qi::rule<std::string::const_iterator, void()> parameter_hunits;
                    qi::rule<std::string::const_iterator, void()> parameter_units;
                    qi::rule<std::string::const_iterator, void()> parameter_temp;
                    qi::rule<std::string::const_iterator, void()> parameter_mass;
                    qi::rule<std::string::const_iterator, void()> parameter_viscosity;
                    qi::rule<std::string::const_iterator, void()> parameter_bf;
                    qi::rule<std::string::const_iterator, void()> parameter;
                    qi::rule<std::string::const_iterator, void()> start;
            };
            // End class BodParserGrammar
    };
    // End class BodParser
}

#endif
