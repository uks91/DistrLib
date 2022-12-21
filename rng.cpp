#include "rng.h"
#include <chrono>
#include <random>

/*unsigned long RNG::mt[RNG_P];
int RNG::mti=RNG_P + 1;


RNG::RNG(unsigned long s)
{
    init_genrand(s);
}

void RNG::init_genrand (unsigned long s)
{
    if (s == 0)
        s = QDateTime::currentMSecsSinceEpoch();


    mt[0] = s & 0xffffffffUL;
    for (mti = 1; mti < RNG_P; mti++) {
        mt[mti] = (1664525UL * mt[mti-1] + 1UL);
        mt[mti]&= 0xffffffffUL;
    }
}

unsigned long RNG::genrand ()
{
    unsigned long y;
    static unsigned long mag01[2] = {0x0UL, RNG_MATRIX_A};

    if (mti >= RNG_P) {
        int kk;

        if (mti == RNG_P+1)
            init_genrand(5489UL);

        for (kk=0; kk<RNG_P-RNG_Q; kk++) {
            y = (mt[kk]&RNG_UPPER_MASK) | (mt[kk+1] & RNG_LOWER_MASK);
            mt[kk] = mt[kk + RNG_Q]^(y>>1)^mag01[y&0x1UL];
        }

        for (; kk<RNG_P-1; kk++) {
            y = (mt[kk]&RNG_UPPER_MASK)|(mt[kk+1]&RNG_LOWER_MASK);
            mt[kk] = mt[kk+(RNG_Q - RNG_P)]^(y>>1)^mag01[y&0x1UL];
        }

        y = (mt[RNG_P-1]&RNG_UPPER_MASK) | (mt[0] & RNG_LOWER_MASK);
        mt[RNG_P-1] = mt[RNG_Q - 1]^(y>>1)^mag01[y&0x1UL];
        mti = 0;
    }

    y = mt[mti++];
    y ^=(y >> 11);
    y ^=(y << 7)  & 0x9d2c5680UL;
    y ^=(y << 15) & 0xefc60000UL;
    y ^=(y >> 18);

    return y;
}

double RNG::uniform(double from, double till)
{
    numGenerations++;
    return from + genrand()*1.0/RNG_MAX_RAND*1.0*(till - from);
}

long RNG::genrand_31()
{
    return (long)(genrand() >> 1);
}*/

int RNG::numGenerations = 0;

class RNG::Impl
{
public:

    double uniform (double from, double till) const
    {
        auto distr = std::uniform_real_distribution<double> (from, till);
        return distr (engine);
    }
//private:
    static std::mt19937 engine;
    static long long timeSinceStart ()
    {
        using namespace std::chrono;
        auto now = steady_clock::now();
        return duration_cast<microseconds> (now.time_since_epoch()).count();
    }
};

std::mt19937 RNG::Impl::engine {std::mt19937(RNG::Impl::timeSinceStart())};

RNG::RNG()
    : impl{new Impl()}
{

}


RNG::~RNG()
{
    delete impl;
}

double RNG::uniform(double from, double till) const
{
    ++numGenerations;
    return impl->uniform(from, till);
}

void RNG::setSeed(long long seed)
{
    RNG::Impl::engine = std::mt19937(seed);
}
