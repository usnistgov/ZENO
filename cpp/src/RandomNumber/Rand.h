#ifndef RAND_H_
#define RAND_H_

#include <random>

class Rand
{
public:
  Rand(int streamNum, int numStreams, int seed);
  ~Rand();

  double getRandIn01();
  double getRandInRange(double min, double max);
  
  // offset
  int k;
  // number of streams
  int m;
  // seed
  int s;
  
  // generator
  std::mt19937_64 gen;
  // distribution
  std::uniform_real_distribution<double> dist;
};

#endif

