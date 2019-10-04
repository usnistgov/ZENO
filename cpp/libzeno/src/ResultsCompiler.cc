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
// Created: Wed Apr 22 11:11:48 2015 EDT
//
// ================================================================

#ifdef USE_MPI
#include <mpi.h>
#endif

#include <cassert>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <thread>

#include "ResultsCompiler.h"

// ================================================================

using namespace zeno;

ResultsCompiler::
ResultsCompiler(ParametersResults const & parameters) 
  : parameters(&parameters) {

}

ResultsCompiler::
~ResultsCompiler() {

}

/// Derives physical quantities from the Walk-on-Spheres and Interior Sampling
/// results.  Either of the results may be NULL, in which case only some of the
/// physical quantities will be computed.  Computation of Form Factors is 
/// optional, due to speed.
///
void 
ResultsCompiler::compile(ResultsZeno const * resultsZeno,
			 ResultsInterior const * resultsInterior,
			 bool computeForm,
			 Results * results) {

  if ((resultsZeno != NULL) &&
      (resultsZeno->getNumHits() == 0.)) {

    std::cerr << "*** Warning ***" << std::endl
	      << "Number of walker hits is zero.  "
	      << "Corresponding results will not be computed.  "
	      << "Try increasing the number of walks?" << std::endl
	      << std::endl;

    resultsZeno = NULL;
  }

  if ((resultsInterior != NULL) &&
      (resultsInterior->getNumHits() == 0.)) {

    std::cerr << "*** Warning ***" << std::endl
	      << "Number of interior sampler hits is zero.  "
	      << "Corresponding results will not be computed.  "
	      << "Try increasing the number of interior samples?" << std::endl
	      << std::endl;

    resultsInterior = NULL;
  }

  if (resultsZeno != NULL) {

    double numWalks = resultsZeno->getNumWalks();

    results->numWalks.set("Number of walks performed",
			  "num_walks_performed",
			  (long long)numWalks,
			  "1");

    Uncertain<double> numZenoHits = resultsZeno->getNumHits();

    Vector3<Uncertain<double> > KPlus  = resultsZeno->getKPlus();
    Vector3<Uncertain<double> > KMinus = resultsZeno->getKMinus();

    Matrix3x3<Uncertain<double> > VPlus  = resultsZeno->getVPlus();
    Matrix3x3<Uncertain<double> > VMinus = resultsZeno->getVMinus();

    Uncertain<double>             t = numZenoHits/numWalks;
    Vector3<Uncertain<double> >   u = (KPlus - KMinus)/numWalks;
    Matrix3x3<Uncertain<double> > v = (VPlus + VMinus)/numWalks;
    Matrix3x3<Uncertain<double> > w = (VPlus - VMinus)/numWalks;

    results->t.set("t", "t", t, "1");
    results->u.set("u", "u", u, "1");
    results->v.set("v", "v", v,
		   Units::getName(parameters->getLengthScaleUnit()));
    results->w.set("w", "w", w,
		   Units::getName(parameters->getLengthScaleUnit()));

    results->capacitance =
      computeCapacitance(t,
			 resultsZeno->getBoundingSphere().getRadius());

    results->polarizabilityTensor =
      computePolarizability(t, u, v, w, 
			    resultsZeno->getBoundingSphere().getRadius());

    results->meanPolarizability =
      computeMeanPolarizability(results->polarizabilityTensor.value);

    Vector3<Uncertain<double> > polarizabilityEigenvalues;
    results->polarizabilityTensor.value.getEigenValues
      (polarizabilityEigenvalues);

    results->polarizabilityEigenvalues.set
      ("Eigenvalues of electric polarizability tensor",
       "electric_polarizability_eigenvalues",
       polarizabilityEigenvalues,
       Units::getName(parameters->getLengthScaleUnit()) + "^3");

    results->hydrodynamicRadius =
      computeHydrodynamicRadius(results->capacitance.value);

    results->q_eta =
      computePadeApproximant(results->polarizabilityTensor.value);

    results->viscometricRadius =
      computeViscometricRadius(results->meanPolarizability.value,
			       results->q_eta.value);

    if (parameters->getMassWasSet()) {

      results->intrinsicViscosityConventional = 
	computeIntrinsicViscosityConventional(results->q_eta.value,
					      results->meanPolarizability.value,
					      parameters->getMassNumber());

      results->intrinsicViscosityConventionalComputed = true;
    }

    if ((parameters->getLengthScaleUnit() != Units::Length::L) &&
	parameters->getSolventViscosityWasSet()) {

      results->frictionCoefficient = 
	computeFrictionCoefficient(parameters->getSolventViscosityNumber(),
				   results->hydrodynamicRadius.value);

      results->frictionCoefficientComputed = true;
    }

    if ((parameters->getLengthScaleUnit() != Units::Length::L) &&
	parameters->getSolventViscosityWasSet() &&
	parameters->getTemperatureWasSet()) {

      results->diffusionCoefficient = 
	computeDiffusionCoefficient(parameters->getTemperatureNumber(),
				    parameters->getSolventViscosityNumber(),
				    results->hydrodynamicRadius.value);

      results->diffusionCoefficientComputed = true;
    }

    if ((parameters->getLengthScaleUnit() != Units::Length::L) &&
	parameters->getSolventViscosityWasSet() &&
	parameters->getBuoyancyFactorWasSet() &&
	parameters->getMassWasSet()) {

      results->sedimentationCoefficient = 
	computeSedimentationCoefficient(parameters->getMassNumber(),
					parameters->getBuoyancyFactor(),
					parameters->getSolventViscosityNumber(),
					results->hydrodynamicRadius.value);

      results->sedimentationCoefficientComputed = true;
    }

    results->resultsZenoCompiled = true;
  }

  if (resultsInterior != NULL) {

    double boundingSphereVolume =
      resultsInterior->getBoundingSphere().getVolume();

    double numInteriorSamples = resultsInterior->getNumSamples();

    results->numInteriorSamples.set("Number of interior samples performed",
				    "num_interior_samples_performed",
				    (long long)numInteriorSamples,
				    "1");
    
    Uncertain<double> numInteriorHits = resultsInterior->getNumHits();

    results->numInteriorHits.set("Interior hits",
				 "interior_hits",
				 numInteriorHits,
				 "1");

    results->volume = computeVolume(numInteriorHits,
				    numInteriorSamples,
				    boundingSphereVolume);

    results->capacitanceOfASphere =
      computeCapacitanceOfASphere(results->volume.value);

    Matrix3x3<Uncertain<double> > hitPointsSqrSum = 
      resultsInterior->getHitPointsSqrSum();

    Vector3<Uncertain<double> > hitPointsSum = 
      resultsInterior->getHitPointsSum();

    results->gyrationTensor = computeGyrationTensor(hitPointsSqrSum, 
						    hitPointsSum,
						    numInteriorHits);

    Vector3<Uncertain<double> > gyrationEigenvalues;
    results->gyrationTensor.value.getEigenValues(gyrationEigenvalues);

    results->gyrationEigenvalues.set
      ("Eigenvalues of gyration tensor",
       "gyration_eigenvalues",
       gyrationEigenvalues,
       Units::getName(parameters->getLengthScaleUnit()) + "^2");
    
    results->resultsInteriorCompiled = true;
  }

  if (resultsZeno != NULL &&
      resultsInterior != NULL) {

    results->intrinsicConductivity =
      computeIntrinsicConductivity(results->meanPolarizability.value,
				   results->volume.value);

    results->intrinsicViscosity = 
      computeIntrinsicViscosity(results->q_eta.value,
				results->intrinsicConductivity.value);
  }
}

Result<Uncertain<double> > 
ResultsCompiler::
computeCapacitance(Uncertain<double> const & t, 
		   double boundingSphereRadius) const {

  Uncertain<double> capacitance = t * boundingSphereRadius;

  const double l = parameters->getLengthScaleNumber();

  capacitance *= l;

  std::string unit = Units::getName(parameters->getLengthScaleUnit());
  
  Result<Uncertain<double> >
    result("Capacitance",
	   "capacitance",
	   capacitance,
	   unit);
  
  return result;
}

Result<Matrix3x3<Uncertain<double> > >
ResultsCompiler::
computePolarizability(Uncertain<double> const & t,
		      Vector3<Uncertain<double> > const & u,
		      Matrix3x3<Uncertain<double> > const & v,
		      Matrix3x3<Uncertain<double> > const & w,
		      double boundingSphereRadius) const {

  Matrix3x3<Uncertain<double> > polarizabilityTensor;

  for (int row = 0; row < 3; row++) {
    for (int col = 0; col < 3; col++) {
      Uncertain<double> element = 
	12 * M_PI * std::pow(boundingSphereRadius, 2)*
	(w.get(row, col) - u.get(row)*v.get(row, col)/t);

      polarizabilityTensor.set(row, col, element);
    }
  }

  polarizabilityTensor.symmetrize();

  const double l = parameters->getLengthScaleNumber();

  polarizabilityTensor *= std::pow(l, 3);

  Result<Matrix3x3<Uncertain<double> > >
    result("Electric polarizability tensor",
	   "electric_polarizability",
	   polarizabilityTensor,
	   Units::getName(parameters->getLengthScaleUnit()) + "^3");

  return result;
}

Result<Uncertain<double> >
ResultsCompiler::
computeMeanPolarizability(Matrix3x3<Uncertain<double> > const & 
			  polarizabilityTensor) 
  const {
  
  Uncertain<double> meanPolarizability = 
    (polarizabilityTensor.get(0, 0) +
     polarizabilityTensor.get(1, 1) +
     polarizabilityTensor.get(2, 2)) / 3.;

  Result<Uncertain<double> >
    result("Mean electric polarizability",
	   "mean_electric_polarizability",
	   meanPolarizability,
	   Units::getName(parameters->getLengthScaleUnit()) + "^3");

  return result;
}

Result<Uncertain<double> > 
ResultsCompiler::
computePadeApproximant(Matrix3x3<Uncertain<double> > const & 
		       polarizabilityTensor)
  const {

  Uncertain<double> alpha1 = polarizabilityTensor.get(0, 0);
  Uncertain<double> alpha2 = polarizabilityTensor.get(1, 1);
  Uncertain<double> alpha3 = polarizabilityTensor.get(2, 2);

  //sort
  if (alpha2 < alpha1) std::swap(alpha1, alpha2);
  if (alpha3 < alpha1) std::swap(alpha1, alpha3);
  if (alpha3 < alpha2) std::swap(alpha2, alpha3);

  Uncertain<double> alpha2_alpha1 = alpha2/alpha1;
  Uncertain<double> alpha3_alpha2 = alpha3/alpha2;

  if ((alpha2_alpha1 < 0.) ||
      (alpha3_alpha2 < 0.)) {

    std::cerr << "*** Warning ***" << std::endl
	      << "Could not compute Prefactor for computing intrinsic " 
	      << "viscosity, Viscometric radius, or Intrinsic viscosity.  "
	      << "This is likely due to Electric polarizability tensor "
	      << "standard deviations being too high.  "
	      << "Try increasing the number of walks?" << std::endl
	      << std::endl;
  }

  Uncertain<double> x1 = log(alpha2_alpha1);
  Uncertain<double> x2 = log(alpha3_alpha2);

  Uncertain<double> q_eta = computePadeApproximant(x1, x2);

  q_eta = Uncertain<double>(q_eta.getMean(),
			    std::pow(q_eta.getMean() * 0.015, 2));

  Result<Uncertain<double> >
    result("Prefactor for computing intrinsic viscosity",
	   "intrinsic_viscosity_prefactor",
	   q_eta,
	   "1");

  return result;
}

Uncertain<double> 
ResultsCompiler::
computePadeApproximant(Uncertain<double> x1, Uncertain<double> x2) const {
  const double delta_i[4] = {4.8,    0.66,  -1.247,  0.787};
  const double k_i[4]     = {0,      1.04,   2.012,  2.315};
  const double b_i[4]     = {0.68,  -7.399,  1.048,  0.136};
  const double t_i[4]     = {0,      1.063,  0.895,  4.993};
  const double B_i[4]     = {1.925, -8.611,  1.652, -0.120};
  const double q_i[4]     = {0,      1.344,  2.029,  1.075};
  const double c_i[4]     = {13.43,  16.17,  0.51,  -5.86};
  const double r_i[4]     = {0,      0.489,  0.879,  2.447};
  const double A_i[4]     = {16.23, -15.92,  14.83, -3.74};
  const double v_i[4]     = {0,      0.462,  1.989,  4.60};
  const double m_i[4]     = {2.786,  0.293, -0.11,   0.012};
  const double u_i[4]     = {0,      0.556,  2.034,  3.024};

  Uncertain<double> delta = 0;
  Uncertain<double> b     = 0;
  Uncertain<double> B     = 0;
  Uncertain<double> c     = 0;
  Uncertain<double> A     = 0;
  Uncertain<double> m     = 0;

  for (int i = 0; i < 4; i++) {
    delta += delta_i[i] * exp(-k_i[i] * x1);
    b     += b_i[i]     * exp(-t_i[i] * x1);
    B     += B_i[i]     * exp(-q_i[i] * x1);
    c     += c_i[i]     * exp(-r_i[i] * x1);
    A     += A_i[i]     * exp(-v_i[i] * x1);
    m     += m_i[i]     * exp(-u_i[i] * x1);
  }

  Uncertain<double> q_eta = 
    (delta*A + c*x2 + b*pow(x2, 2.) + 4.*pow(x2, m)) /
    (6.*A + 6.*c*x2/delta + B*pow(x2, 2.) + 5.*pow(x2, m));

  return q_eta;
}

Result<Uncertain<double> >
ResultsCompiler::
computeVolume(Uncertain<double> const & numInteriorHits, 
	      double numInteriorSamples,
	      double boundingSphereVolume) const {

  Uncertain<double> volume = 
    boundingSphereVolume * numInteriorHits / numInteriorSamples;

  const double l = parameters->getLengthScaleNumber();

  volume *= std::pow(l, 3);

  Result<Uncertain<double> >
    result("Volume",
	   "volume",
	   volume,
	   Units::getName(parameters->getLengthScaleUnit()) + "^3");
  
  return result;
}

Result<Uncertain<double> >
ResultsCompiler::
computeIntrinsicConductivity(Uncertain<double> const & meanPolarizability,
			     Uncertain<double> const & volume) const {

  Uncertain<double> intrinsicConductivity = meanPolarizability / volume;

  Result<Uncertain<double> >
    result("Intrinsic conductivity",
	   "intrinsic_conductivity",
	   intrinsicConductivity,
	   "1");
		
  
  return result;
}

Result<Uncertain<double> > 
ResultsCompiler::
computeCapacitanceOfASphere(Uncertain<double> const & volume) const {

  Uncertain<double> capacitanceOfASphere = pow(3.*volume/(4.*M_PI), 1./3.);

  Result<Uncertain<double> >
    result("Capacitance of a sphere of the same volume",
	   "capacitance_sphere_same_volume",
	   capacitanceOfASphere,
	   Units::getName(parameters->getLengthScaleUnit()));
  
  return result;
}

Result<Uncertain<double> > 
ResultsCompiler::
computeHydrodynamicRadius(Uncertain<double> const & capacitance) const {

  Uncertain<double> q_Rh(1, std::pow(0.01, 2));

  Uncertain<double> hydrodynamicRadius = q_Rh * capacitance;

  Result<Uncertain<double> >
    result("Hydrodynamic radius",
	   "hydrodynamic_radius",
	   hydrodynamicRadius,
	   Units::getName(parameters->getLengthScaleUnit()));

  return result;
}

Result<Uncertain<double> > 
ResultsCompiler::
computeViscometricRadius(Uncertain<double> const & meanPolarizability,
			 Uncertain<double> const & padeApproximant) const {

  Uncertain<double> alpha = meanPolarizability;
  Uncertain<double> q_eta = padeApproximant;

  Uncertain<double> viscometricRadius = pow((3.*q_eta*alpha)/(10.*M_PI), 1./3.);

  Result<Uncertain<double> >
    result("Viscometric radius",
	   "viscometric_radius",
	   viscometricRadius,
	   Units::getName(parameters->getLengthScaleUnit()));
  
  return result;
}

Result<Uncertain<double> > 
ResultsCompiler::
computeIntrinsicViscosity(Uncertain<double> const & padeApproximant,
			  Uncertain<double> const & intrinsicConductivity) 
  const {

  Uncertain<double> intrinsicViscosity = 
    padeApproximant * intrinsicConductivity;

  Result<Uncertain<double> >
    result("Intrinsic viscosity",
	   "intrinsic_viscosity",
	   intrinsicViscosity,
	   "1");
  
  return result;
}

Result<Uncertain<double> >
ResultsCompiler::
computeIntrinsicViscosityConventional
(Uncertain<double> const & padeApproximant,
 Uncertain<double> const & meanPolarizability,
 double mass) const {
  
  Uncertain<double> intrinsicViscosityConventional =
    padeApproximant * meanPolarizability / mass;

  std::string units =
    Units::getName(parameters->getLengthScaleUnit()) +
    "^3/" +
    Units::getName(parameters->getMassUnit());
  
  if (parameters->getLengthScaleUnit() != Units::Length::L) {
    const Uncertain<double> a_l = 
      Units::getFactor(parameters->getLengthScaleUnit(), 
		       Units::Length::cm);

    const Uncertain<double> a_m = 
      Units::getFactor(parameters->getMassUnit(), 
		       Units::Mass::g);

    intrinsicViscosityConventional *= pow(a_l, 3.) * a_m;

    units = "cm^3/g";
  }

  Result<Uncertain<double> >
    result("Intrinsic viscosity with mass units",
	   "intrinsic_viscosity_mass_units",
	   intrinsicViscosityConventional,
	   units);

  return result;
}

Result<Uncertain<double> >
ResultsCompiler::
computeFrictionCoefficient(double solventViscosity,
			   Uncertain<double> const & hydrodynamicRadius) const {

  Uncertain<double> frictionCoefficient = 
    6 * M_PI * solventViscosity * hydrodynamicRadius;

  const Uncertain<double> a_l =
    Units::getFactor(parameters->getLengthScaleUnit(),
                     Units::Length::cm);

  const Uncertain<double> a_eta =
    Units::getFactor(parameters->getSolventViscosityUnit(),
		     Units::Viscosity::cp);

  frictionCoefficient *= std::pow(10, -2) * a_l * a_eta;

  Result<Uncertain<double> >
    result("Friction coefficient",
	   "friction_coefficient",
	   frictionCoefficient,
	   "d.s/cm");
  
  return result;
}

Result<Uncertain<double> >
ResultsCompiler::
computeDiffusionCoefficient(double temperature,
			    double solventViscosity,
			    Uncertain<double> const & hydrodynamicRadius) 
  const {

  const Uncertain<double> a_T =
    Units::getOffset(parameters->getTemperatureUnit(),
		     Units::Temperature::K);

  const Uncertain<double> k_B = Units::kB();

  Uncertain<double> diffusionCoefficient = 
    k_B * (temperature + a_T) / 
    (6 * M_PI * solventViscosity * hydrodynamicRadius);

  const Uncertain<double> a_l =
    Units::getFactor(parameters->getLengthScaleUnit(),
                     Units::Length::cm);

  const Uncertain<double> a_eta =
    Units::getFactor(parameters->getSolventViscosityUnit(),
		     Units::Viscosity::cp);

  diffusionCoefficient *= 
    std::pow(10, 9) * pow(a_l, -1.) * pow(a_eta, -1.);

  Result<Uncertain<double> >
    result("Diffusion coefficient",
	   "diffusion_coefficient",
	   diffusionCoefficient,
	   "cm^2/s");
  
  return result;
}

Result<Uncertain<double> >
ResultsCompiler::
computeSedimentationCoefficient(double mass,
				double buoyancyFactor,
				double solventViscosity,
				Uncertain<double> const & hydrodynamicRadius) 
  const {

  Uncertain<double> sedimentationCoefficient = 
    mass * buoyancyFactor / 
    (6 * M_PI * solventViscosity * hydrodynamicRadius);

  const Uncertain<double> a_l =
    Units::getFactor(parameters->getLengthScaleUnit(),
                     Units::Length::cm);

  const Uncertain<double> a_eta =
    Units::getFactor(parameters->getSolventViscosityUnit(),
		     Units::Viscosity::cp);

  const Uncertain<double> a_m =
    Units::getFactor(parameters->getMassUnit(),
		     Units::Mass::g);

  sedimentationCoefficient *= 
    std::pow(10, 15) *
    pow(a_l, -1.) * pow(a_eta, -1.) * pow(a_m, -1.);

  Result<Uncertain<double> >
    result("Sedimentation coefficient",
	   "sedimentation_coefficient",
	   sedimentationCoefficient,
	   "Sved");
  
  return result;
}

Result<Matrix3x3<Uncertain<double> > >
ResultsCompiler::
computeGyrationTensor(Matrix3x3<Uncertain<double> > const & hitPointsSqrSum,
		      Vector3<Uncertain<double> > const & hitPointsSum,
		      Uncertain<double> const & numInteriorHits) const {

  Matrix3x3<Uncertain<double> > hitPointsSumSqr;

  for (int row = 0; row < 3; ++row) {
    for (int col = 0; col < 3; ++col) {
      Uncertain<double> element = 
	hitPointsSum.get(row) * hitPointsSum.get(col);

      hitPointsSumSqr.set(row, col, element);
    }
  }

  Matrix3x3<Uncertain<double> > gyrationTensor =
    hitPointsSqrSum / numInteriorHits - 
    hitPointsSumSqr / pow(numInteriorHits, 2.);

  const double l = parameters->getLengthScaleNumber();

  gyrationTensor *= std::pow(l, 2);

  Result<Matrix3x3<Uncertain<double> > >
    result("Gyration tensor",
	   "gyration",
	   gyrationTensor,
	   Units::getName(parameters->getLengthScaleUnit()) + "^2");
  
  return result;
}

/// Inverts the formula: 
/// index = (j*j + j)/2 + i
///
void 
ResultsCompiler::
indexToIJ(BigUInt index, BigUInt * i, BigUInt * j) const {

  (*j) = round(std::sqrt(2*index + 1)) - 1;
  (*i) = index - ((*j)*(*j) + (*j))/2;

  ++(*j);
}

void
ResultsCompiler::
printScalar(std::string const & prettyName,
	    std::string const & csvName,
	    std::string const & units,
	    Uncertain<double> const & property,
	    std::ofstream * csvOutputFile) const {

  std::cout << std::scientific
	    << prettyName << " (" << units << "): " << property << std::endl
	    << std::endl;
  
  *csvOutputFile << std::scientific
		 << csvName << ",units," << units
		 << std::endl
		 << csvName << ",mean," << property.getMean()
		 << std::endl
		 << csvName << ",std_dev," << property.getStdDev()
		 << std::endl;
}

void
ResultsCompiler::
printVector3(std::string const & prettyName,
	     std::string const & csvName,
	     std::string const & units,
	     Vector3<Uncertain<double> > const & property,
	     std::ofstream * csvOutputFile) const {
  
  std::cout << std::scientific
	    << prettyName << " (" << units << "): " << std::endl
	    << property << std::endl
	    << std::endl;
  
  for (int i = 0; i < 3; ++i) {
    *csvOutputFile << std::scientific
		   << csvName << "[" << i << "],"
		   << "units,"
		   << units
		   << std::endl
		   << csvName << "[" << i << "],"
		   << "mean,"
		   << property.get(i).getMean()
		   << std::endl
		   << csvName << "[" << i << "],"
		   << "std_dev,"
		   << property.get(i).getStdDev()
		   << std::endl;
  }
}

void
ResultsCompiler::
printMatrix3x3(std::string const & prettyName,
	       std::string const & csvName,
	       std::string const & units,
	       Matrix3x3<Uncertain<double> > const & property,
	       std::ofstream * csvOutputFile) const {

  std::cout << std::scientific
	    << prettyName << " (" << units << "): " << std::endl
	    << property << std::endl
	    << std::endl;
  
  for (int row = 0; row < 3; ++row) {
    for (int col = 0; col < 3; ++col) {
      *csvOutputFile << std::scientific
		     << csvName << "[" << row << "][" << col << "],"
		     << "units,"
		     << units
		     << std::endl
		     << csvName << "[" << row << "][" << col << "],"
		     << "mean,"
		     << property.get(row, col).getMean()
		     << std::endl
		     << csvName << "[" << row << "][" << col << "],"
		     << "std_dev,"
		     << property.get(row, col).getStdDev()
		     << std::endl;
    }
  }
}
