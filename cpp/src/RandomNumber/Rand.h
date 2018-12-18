#ifndef RAND_H_
#define RAND_H_

class Rand
{
public:
  Rand(int streamNum, int numStreams, int seed);
  ~Rand();

  double getRandIn01();
  double getRandInRange(double min, double max);
};

#endif

