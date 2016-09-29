// ================================================================
// 
// Disclaimer:  IMPORTANT:  This software was developed at the
// National Institute of Standards and Technology by employees of the
// Federal Government in the course of their official duties.
// Pursuant to title 17 Section 105 of the United States Code this
// software is not subject to copyright protection and is in the
// public domain.  This is an experimental system.  NIST assumes no
// responsibility whatsoever for its use by other parties, and makes
// no guarantees, expressed or implied, about its quality,
// reliability, or any other characteristic.  We would appreciate
// acknowledgement if the software is used.  This software can be
// redistributed and/or modified freely provided that any derivative
// works bear some notice that they are derived from it, and any
// modified versions bear some notice that they have been modified.
// 
// ================================================================

// ================================================================
// 
// Authors: Derek Juba <derek.juba@nist.gov>
// Date:    Wed Apr 22 11:11:48 2015 EDT
//
// Time-stamp: <2016-09-28 12:07:33 dcj>
//
// ================================================================

#ifdef USE_MPI
#include <mpi.h>
#endif

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <thread>

#include "ResultsCompiler.h"

// ================================================================

ResultsCompiler::
ResultsCompiler(Parameters const & parameters) 
  : parameters(&parameters),
    boundingSphereRadius(),
    boundingSphereCenter(),
    t(),
    u(),
    v(),
    w(),
    capacitance(),
    polarizabilityTensor(),
    meanPolarizability(),
    polarizabilityEigenvalues(),
    volume(),
    intrinsicConductivity(),
    capacitanceOfASphere(),
    hydrodynamicRadius(),
    q_eta(),
    viscometricRadius(),
    intrinsicViscosity(),
    intrinsicViscosityConventional(),
    frictionCoefficient(),
    diffusionCoefficient(),
    sedimentationCoefficient(),
    gyrationTensor(),
    gyrationEigenvalues(),
    formFactors(),
    intrinsicViscosityConventionalComputed(false),
    frictionCoefficientComputed(false),
    diffusionCoefficientComputed(false),
    sedimentationCoefficientComputed(false),
    resultsZenoCompiled(false), 
    resultsInteriorCompiled(false),
    formResultsCompiled(false) {

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
			 Sphere<double> const & boundingSphere,
			 bool computeForm) {

  boundingSphereRadius = boundingSphere.getRadius();
  boundingSphereCenter = boundingSphere.getCenter();

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

    Uncertain<double> numZenoHits = resultsZeno->getNumHits();

    Vector3<Uncertain<double> > KPlus  = resultsZeno->getKPlus();
    Vector3<Uncertain<double> > KMinus = resultsZeno->getKMinus();

    Matrix3x3<Uncertain<double> > VPlus  = resultsZeno->getVPlus();
    Matrix3x3<Uncertain<double> > VMinus = resultsZeno->getVMinus();

    t = numZenoHits/numWalks;
    u = (KPlus - KMinus)/numWalks;
    v = (VPlus + VMinus)/numWalks;
    w = (VPlus - VMinus)/numWalks;

    capacitance = computeCapacitance(t, boundingSphereRadius);

    polarizabilityTensor = computePolarizability(t, u, v, w, 
						 boundingSphereRadius);

    meanPolarizability = computeMeanPolarizability(polarizabilityTensor);

    polarizabilityTensor.getEigenValues(polarizabilityEigenvalues);

    hydrodynamicRadius = computeHydrodynamicRadius(capacitance);

    q_eta = computePadeApproximant(polarizabilityTensor);

    viscometricRadius = computeViscometricRadius(meanPolarizability,
						 q_eta);

    if (parameters->getMassWasSet()) {

      intrinsicViscosityConventional = 
	computeIntrinsicViscosityConventional(q_eta,
					      meanPolarizability,
					      parameters->getMassNumber());

      intrinsicViscosityConventionalComputed = true;
    }

    if ((parameters->getLengthScaleUnit() != Units::Length::L) &&
	parameters->getSolventViscosityWasSet()) {

      frictionCoefficient = 
	computeFrictionCoefficient(parameters->getSolventViscosityNumber(),
				   hydrodynamicRadius);

      frictionCoefficientComputed = true;
    }

    if ((parameters->getLengthScaleUnit() != Units::Length::L) &&
	parameters->getSolventViscosityWasSet() &&
	parameters->getTemperatureWasSet()) {

      diffusionCoefficient = 
	computeDiffusionCoefficient(parameters->getTemperatureNumber(),
				    parameters->getSolventViscosityNumber(),
				    hydrodynamicRadius);

      diffusionCoefficientComputed = true;
    }

    if ((parameters->getLengthScaleUnit() != Units::Length::L) &&
	parameters->getSolventViscosityWasSet() &&
	parameters->getBuoyancyFactorWasSet() &&
	parameters->getMassWasSet()) {

      sedimentationCoefficient = 
	computeSedimentationCoefficient(parameters->getMassNumber(),
					parameters->getBuoyancyFactor(),
					parameters->getSolventViscosityNumber(),
					hydrodynamicRadius);

      sedimentationCoefficientComputed = true;
    }

    resultsZenoCompiled = true;
  }

  if (resultsInterior != NULL) {

    double boundingSphereVolume = boundingSphere.getVolume();

    double numInteriorSamples = resultsInterior->getNumSamples();

    numInteriorHits = resultsInterior->getNumHits();

    volume = computeVolume(numInteriorHits, numInteriorSamples,
			   boundingSphereVolume);

    capacitanceOfASphere = computeCapacitanceOfASphere(volume);

    Matrix3x3<Uncertain<double> > hitPointsSqrSum = 
      resultsInterior->getHitPointsSqrSum();

    Vector3<Uncertain<double> > hitPointsSum = 
      resultsInterior->getHitPointsSum();

    gyrationTensor = computeGyrationTensor(hitPointsSqrSum, 
					   hitPointsSum,
					   numInteriorHits);

    gyrationTensor.getEigenValues(gyrationEigenvalues);

    resultsInteriorCompiled = true;
  }

  if (resultsZeno != NULL &&
      resultsInterior != NULL) {

    intrinsicConductivity = computeIntrinsicConductivity(meanPolarizability,
							 volume);

    intrinsicViscosity = 
      computeIntrinsicViscosity(q_eta, intrinsicConductivity);
  }

  if (computeForm) {

    assert(resultsInterior != NULL);

    formFactorQs = computeFormFactorQs(boundingSphereRadius);

    formFactors = computeFormFactors(*(resultsInterior->getPoints()),
					       formFactorQs);

    formResultsCompiled = true;
  }
}

Uncertain<double> 
ResultsCompiler::
computeCapacitance(Uncertain<double> const & t, 
		   double boundingSphereRadius) const {

  Uncertain<double> capacitance = t * boundingSphereRadius;

  const double l = parameters->getLengthScaleNumber();

  capacitance *= l;

  return capacitance;
}

Matrix3x3<Uncertain<double> >
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
	12*M_PI*pow(boundingSphereRadius, 2)*
	(w.get(row, col) - u.get(row)*v.get(row, col)/t);

      polarizabilityTensor.set(row, col, element);
    }
  }

  polarizabilityTensor.symmetrize();

  const double l = parameters->getLengthScaleNumber();

  polarizabilityTensor *= pow(l, 3);

  return polarizabilityTensor;
}

Uncertain<double>
ResultsCompiler::
computeMeanPolarizability(Matrix3x3<Uncertain<double> > const & 
			  polarizabilityTensor) 
  const {
  
  Uncertain<double> meanPolarizability = 
    (polarizabilityTensor.get(0, 0) +
     polarizabilityTensor.get(1, 1) +
     polarizabilityTensor.get(2, 2)) / 3.;

  return meanPolarizability;
}

Uncertain<double> 
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
	      << "This is likely due to Electrostatic polarizability tensor "
	      << "standard deviations being too high.  "
	      << "Try increasing the number of walks?" << std::endl
	      << std::endl;
  }

  Uncertain<double> x1 = log(alpha2_alpha1);
  Uncertain<double> x2 = log(alpha3_alpha2);

  Uncertain<double> q_eta = computePadeApproximant(x1, x2);

  q_eta = Uncertain<double>(q_eta.getMean(),
			    pow(q_eta.getMean() * 0.015, 2));

  return q_eta;
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

Uncertain<double>
ResultsCompiler::
computeVolume(Uncertain<double> const & numInteriorHits, 
	      double numInteriorSamples,
	      double boundingSphereVolume) const {

  Uncertain<double> volume = 
    boundingSphereVolume * numInteriorHits / numInteriorSamples;

  const double l = parameters->getLengthScaleNumber();

  volume *= pow(l, 3);

  return volume;
}

Uncertain<double>
ResultsCompiler::
computeIntrinsicConductivity(Uncertain<double> const & meanPolarizability,
			     Uncertain<double> const & volume) const {

  Uncertain<double> intrinsicConductivity = meanPolarizability / volume;

  return intrinsicConductivity;
}

Uncertain<double> 
ResultsCompiler::
computeCapacitanceOfASphere(Uncertain<double> const & volume) const {

  Uncertain<double> capacitanceOfASphere = pow(3.*volume/(4.*M_PI), 1./3.);

  return capacitanceOfASphere;
}

Uncertain<double> 
ResultsCompiler::
computeHydrodynamicRadius(Uncertain<double> const & capacitance) const {

  Uncertain<double> q_Rh(1, pow(0.01, 2));

  Uncertain<double> hydrodynamicRadius = q_Rh * capacitance;

  return hydrodynamicRadius;
}

Uncertain<double> 
ResultsCompiler::
computeViscometricRadius(Uncertain<double> const & meanPolarizability,
			 Uncertain<double> const & padeApproximant) const {

  Uncertain<double> alpha = meanPolarizability;
  Uncertain<double> q_eta = padeApproximant;

  Uncertain<double> viscometricRadius = pow((3.*q_eta*alpha)/(10.*M_PI), 1./3.);

  return viscometricRadius;
}

Uncertain<double> 
ResultsCompiler::
computeIntrinsicViscosity(Uncertain<double> const & padeApproximant,
			  Uncertain<double> const & intrinsicConductivity) 
  const {

  Uncertain<double> intrinsicViscosity = 
    padeApproximant * intrinsicConductivity;

  return intrinsicViscosity;
}

Uncertain<double>
ResultsCompiler::
computeIntrinsicViscosityConventional
(Uncertain<double> const & padeApproximant,
 Uncertain<double> const & meanPolarizability,
 double mass) const {

  Uncertain<double> intrinsicViscosityConventional =
    padeApproximant * meanPolarizability / mass;

  if (parameters->getLengthScaleUnit() != Units::Length::L) {
    const Uncertain<double> a_l = 
      Units::getFactor(parameters->getLengthScaleUnit(), 
		       Units::Length::cm);

    const Uncertain<double> a_m = 
      Units::getFactor(parameters->getMassUnit(), 
		       Units::Mass::g);

    intrinsicViscosityConventional *= pow(a_l, 3.) * a_m;	
  }

  return intrinsicViscosityConventional;
}

Uncertain<double>
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

  frictionCoefficient *= pow(10, -2) * a_l * a_eta;

  return frictionCoefficient;
}

Uncertain<double>
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
    pow(10, 9) * pow(a_l, -1.) * pow(a_eta, -1.);

  return diffusionCoefficient;
}

Uncertain<double>
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
    pow(10, 15) * pow(a_l, -1.) * pow(a_eta, -1.) * pow(a_m, -1.);

  return sedimentationCoefficient;
}

Matrix3x3<Uncertain<double> >
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

  gyrationTensor *= pow(l, 2);

  return gyrationTensor;
}

std::array<double, ResultsCompiler::numFormFactors>
ResultsCompiler::
computeFormFactorQs(double boundingSphereRadius) const {

  const double l = parameters->getLengthScaleNumber();

  std::array<double, numFormFactors> qs;

  for (unsigned int factorNum = 0; 
       factorNum < numFormFactors; 
       ++factorNum) {

    qs.at(factorNum) =
      0.01 * pow(10000, (double)factorNum/(numFormFactors - 1)) /
      boundingSphereRadius;

    qs.at(factorNum) *= pow(l, -1);
  }

  return qs;
}

std::array<double, ResultsCompiler::numFormFactors>
ResultsCompiler::
computeFormFactors(std::vector<Vector3<double> > const & interiorPoints,
			std::array<double, 
			ResultsCompiler::numFormFactors> const & 
			formFactorQs) const {

  BigUInt numPairs = 
    (interiorPoints.size() * interiorPoints.size() - 
     interiorPoints.size()) / 2;

  BigUInt numPairsInProcess = numPairs;

  int numThreads = parameters->getNumThreads(); 

  std::vector<std::array<double, numFormFactors> > threadsFormFactors;

  threadsFormFactors.reserve(numThreads);

  std::vector<std::thread> threads;

  threads.reserve(numThreads);

  BigUInt threadStartPairIndex = 0;

  for (int threadNum = 0; threadNum < numThreads; ++threadNum) {
    BigUInt numPairsInThread = numPairsInProcess / numThreads;

    if ((unsigned int)threadNum < numPairsInProcess % numThreads) {
      numPairsInThread ++;
    }

    BigUInt threadEndPairIndex = threadStartPairIndex + numPairsInThread;

    threads.emplace_back(&ResultsCompiler::computeFormFactorsThread,
			 this,
			 threadNum,
			 interiorPoints,
			 formFactorQs,
			 threadStartPairIndex,
			 threadEndPairIndex,
			 &(threadsFormFactors[threadNum]));

    threadStartPairIndex += numPairsInThread;
  }

  for (int threadNum = 0; threadNum < numThreads; ++threadNum) {
    threads.at(threadNum).join();
  }

  std::array<double, numFormFactors> formFactorsReduced;

  formFactorsReduced.fill(0);

  for (int threadNum = 0; 
       threadNum < numThreads; 
       ++threadNum) {

    for (unsigned int factorNum = 0; 
	 factorNum < numFormFactors; 
	 ++factorNum) {

      formFactorsReduced.at(factorNum) += 
	threadsFormFactors[threadNum].at(factorNum);
    }
  }

  for (unsigned int factorNum = 0; 
       factorNum < numFormFactors; 
       ++factorNum) {

    formFactorsReduced.at(factorNum) /= numPairs;
  }

  return formFactorsReduced;
}

void
ResultsCompiler::
computeFormFactorsThread(int threadNum,
			      std::vector<Vector3<double> > const & 
			      interiorPoints,
			      std::array<double, 
			      ResultsCompiler::numFormFactors> const & 
			      formFactorQs,
			      BigUInt startPairIndex,
			      BigUInt endPairIndex,
			      std::array<double, numFormFactors> * 
			      threadFormFactors) const {

  const double l = parameters->getLengthScaleNumber();

  threadFormFactors->fill(0);

  for (BigUInt pairIndex = startPairIndex;
       pairIndex < endPairIndex;
       ++pairIndex) {

    BigUInt i = 0, j = 0;

    indexToIJ(pairIndex, &i, &j);

    double distance = 
      (interiorPoints.at(i) - 
       interiorPoints.at(j)).getMagnitude();

    distance *= l;

    for (unsigned int factorNum = 0; 
	 factorNum < numFormFactors; 
	 ++factorNum) {

      double q = formFactorQs.at(factorNum);

      double divisor = q*distance;

      //use Taylor expansion of sin(x)/x for small x
      if (divisor < 0.001) {
	threadFormFactors->at(factorNum) +=
	  1. - pow(divisor, 2)/6.;
      }
      else {
	threadFormFactors->at(factorNum) += 
	  sin(divisor)/divisor;
      }
    }
  }
}

/// Inverts the formula: 
/// index = (j*j + j)/2 + i
///
void 
ResultsCompiler::
indexToIJ(BigUInt index, BigUInt * i, BigUInt * j) const {

  (*j) = round(sqrt(2*index + 1)) - 1;
  (*i) = index - ((*j)*(*j) + (*j))/2;

  ++(*j);
}

/// Print the computed physical quantities.  Optionally also print the raw
/// hit counts.
///
void 
ResultsCompiler::
print(bool printCounts) const {

  if (resultsZenoCompiled) {

    std::cout << std::scientific
	      << "Capacitance (" 
	      << Units::getName(parameters->getLengthScaleUnit()) 
	      << "): " << capacitance << std::endl
	      << std::endl;

    std::cout << std::scientific
	      << "Electric polarizability tensor ("
	      << Units::getName(parameters->getLengthScaleUnit())
	      << "^3): " << std::endl
	      << polarizabilityTensor 
	      << std::endl;

    std::cout << std::scientific
	      << "Eigenvalues of electric polarizability tensor (" 
	      << Units::getName(parameters->getLengthScaleUnit())
	      << "^3): " << std::endl
	      << polarizabilityEigenvalues << std::endl
	      << std::endl;

    std::cout << std::scientific
	      << "Mean electric polarizability (" 
	      << Units::getName(parameters->getLengthScaleUnit())
	      << "^3): " << meanPolarizability << std::endl
	      << std::endl;

    std::cout << std::scientific
	      << "Hydrodynamic radius ("
	      << Units::getName(parameters->getLengthScaleUnit())
	      << "): " << hydrodynamicRadius << std::endl
	      << std::endl;

    std::cout << std::scientific
	      << "Prefactor for computing intrinsic viscosity: " 
	      << q_eta << std::endl
	      << std::endl;

    std::cout << std::scientific
	      << "Viscometric radius (" 
	      << Units::getName(parameters->getLengthScaleUnit())
	      << "): " << viscometricRadius << std::endl
	      << std::endl;

    if (intrinsicViscosityConventionalComputed) {

      std::cout << "Intrinsic viscosity with mass units (";

      if (parameters->getLengthScaleUnit() == Units::Length::L) {

	std::cout << Units::getName(parameters->getLengthScaleUnit())
		  << "^3/"
		  << Units::getName(parameters->getMassUnit());
      }
      else {

	std::cout << "cm^3/g";
      }

      std::cout << std::scientific
		<< "): " << intrinsicViscosityConventional << std::endl
		<< std::endl;
    }

    if (frictionCoefficientComputed) {

      std::cout << std::scientific
		<< "Friction coefficient (d.s/cm): " 
		<< frictionCoefficient << std::endl
		<< std::endl;
    }

    if (diffusionCoefficientComputed) {

      std::cout << std::scientific
		<< "Diffusion coefficient (cm^2/s): "
		<< diffusionCoefficient << std::endl
		<< std::endl;
    }

    if (sedimentationCoefficientComputed) {

      std::cout << std::scientific
		<< "Sedimentation coefficient (Sved): "
		<< sedimentationCoefficient << std::endl
		<< std::endl;
    }
  }

  if (resultsInteriorCompiled) {

    std::cout << std::scientific
	      << "Volume (" 
	      << Units::getName(parameters->getLengthScaleUnit())
	      << "^3): " << volume << std::endl
	      << std::endl;

    std::cout << std::scientific
	      << "Capacitance of a sphere of the same volume ("
	      << Units::getName(parameters->getLengthScaleUnit())
	      << "): " << capacitanceOfASphere << std::endl
	      << std::endl;

    std::cout << std::scientific
	      << "Gyration tensor ("
	      << Units::getName(parameters->getLengthScaleUnit())
	      << "^2): " << std::endl
	      << gyrationTensor 
	      << std::endl;

    std::cout << std::scientific
	      << "Eigenvalues of gyration tensor (" 
	      << Units::getName(parameters->getLengthScaleUnit())
	      << "^2): " << std::endl
	      << gyrationEigenvalues << std::endl
	      << std::endl;
  }

  if (resultsZenoCompiled && resultsInteriorCompiled) {

    std::cout << std::scientific
	      << "Intrinsic conductivity: " << intrinsicConductivity 
	      << std::endl
	      << std::endl;

    std::cout << std::scientific
	      << "Intrinsic viscosity: " << intrinsicViscosity 
	      << std::endl
	      << std::endl;
  }

  if (formResultsCompiled) {

    std::cout << "Form factor: " << std::endl;

    for (unsigned int factorNum = 0; 
	 factorNum < numFormFactors;
	 ++factorNum) {

      std::cout << std::scientific
		<< formFactorQs.at(factorNum) << " ("
		<< Units::getName(parameters->getLengthScaleUnit())
		<< "^-1): "
		<< formFactors.at(factorNum) << std::endl;
    }

    std::cout << std::endl;
  }

  if (printCounts) {

    std::cout << std::scientific
	      << "Counts:" << std::endl
	      << std::endl
	      << "t: " << t << std::endl
	      << std::endl
	      << "u: " << u << std::endl
	      << std::endl
	      << "v: " << std::endl << v
	      << std::endl
	      << "w: " << std::endl << w
	      << std::endl
	      << "Interior hits: " << numInteriorHits << std::endl
	      << std::endl;
  }
}

Uncertain<double> 
ResultsCompiler::
getCapacitance() const {

  assert(resultsZenoCompiled);

  return capacitance;
}

Uncertain<double> 
ResultsCompiler::
getMeanPolarizability() const {

  assert(resultsZenoCompiled);

  return meanPolarizability;
}

Uncertain<double> 
ResultsCompiler::
getVolume() const {

  assert(resultsInteriorCompiled);

  return volume;
}

// ================================================================

// Local Variables:
// time-stamp-line-limit: 30
// mode: c++
// End:
