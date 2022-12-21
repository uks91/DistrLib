#ifndef RNG_H
#define RNG_H

//#define RNG_P 624
//#define RNG_Q 397
//#define RNG_MATRIX_A 0x9908b0dUL
//#define RNG_UPPER_MASK 0x80000000UL
//#define RNG_LOWER_MASK 0x7fffffffUL
//#define RNG_MAX_RAND 0xffffffffUL

//static unsigned long mt[RNG_P];
//static int mti=RNG_P + 1;

//#include <QDateTime>
//#include "functions.h"


/*class RNG
{
public:
    RNG(unsigned long s=0);

    void init_genrand (unsigned long s=0);
    unsigned long genrand ();
    long genrand_31();
    double uniform(double from, double till);
    double uniform() {return uniform (0,1);}
    static int numGenerations;

private:
    static unsigned long mt[RNG_P];
    static int mti;
};*/

class RNG
{
public:
    RNG();
    ~RNG();
    double uniform (double from=0.0, double till=1.0) const;
    static int numGenerations;
    static void setSeed (long long seed);

private:
    class Impl;
    Impl* impl;
};

#endif // RNG_H
