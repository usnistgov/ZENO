#include "Rand.h"

Rand::Rand(int streamNum, int numStreams, int seed) {
	k = streamNum;
	m = numStreams;
	s = seed;
	
	// initialize generator
	gen = std::mt19937_64(s);
	// offset by k (stagger starting value)
	gen.discard(k);
	// initialize distribution
	dist = std::uniform_real_distribution<double>(0., 1.);
}

Rand::~Rand() {
}

double 
Rand::getRandIn01() {
	// generate
	double value = dist(gen);
	// discard by total number of streams present minus 1
	gen.discard(m-1);
	
	return value;
}

double 
Rand::getRandInRange(double min, double max) {
	// generate
	double value = dist(gen) * (max - min) + min;
	// discard according to total number of streams minus 1
	gen.discard(m-1);
	
	return value;
}

