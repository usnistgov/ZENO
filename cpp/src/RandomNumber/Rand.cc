#include "Rand.h"
#include <random>

// mersenne_twister_engine 64-bit
std::random_device r;
std::seed_seq ran_seed{r(), r(), r(), r(), r(), r(), r(), r()};
std::mt19937_64 gen(ran_seed);

Rand::Rand(int streamNum, int numStreams, int seed) {
	gen.seed(seed);
}

Rand::~Rand() {
}

double 
Rand::getRandIn01() {
	std::uniform_real_distribution<> dist(0.0, 1.0);
	return dist(gen);
}

double 
Rand::getRandInRange(double min, double max) {
	std::uniform_real_distribution<> dist(0.0, 1.0);
	return dist(gen) * (max - min) + min;;
}

