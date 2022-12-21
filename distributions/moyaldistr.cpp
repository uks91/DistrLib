#include "moyaldistr.h"

MoyalDistr::MoyalDistr() : AbstractDistribution()
{
    pars.resize(2);

//    AbstractDistribution::defineParameters();
}

MoyalDistr::MoyalDistr(const array1 &in, const array1 *pars0) : AbstractDistribution(in, pars0)
{
    pars.resize(2);

    AbstractDistribution::defineParameters();
}

qreal MoyalDistr::cdfVerified(qreal x, const array1 &tmpPars)
{
//   qDebug () << tmpPars;
    return erfc ( exp(-(x-tmpPars[MEAN])/(2*tmpPars[SIGMA]))/sqrt(2) );
}

void MoyalDistr::defineParameters(const array1 &d, array1 &pars)
{
    if (pars.size() !=2)
        pars.resize(2);
    pars[MEAN] = getMean(d);
    pars[SIGMA] = getSigma(d, pars[MEAN]);

    pars = estimate(d, pars);
//    qDebug () << getMinOrMax();
}

qreal MoyalDistr::pdf(qreal x, const array1 &tmpPars)
{
    if (tmpPars[SIGMA] <= 0)
        return 0;
    return func( (x-tmpPars[MEAN])/tmpPars[SIGMA] )/tmpPars[SIGMA];
}

qreal MoyalDistr::rnd()
{
    qreal h, z, y;
    do {
        h = rng.uniform(0, 1)*0.912;
        y = M_PI*rng.uniform(0, 1) - M_PI_2;
        z = tan(y);
    } while (h > func(z)/pow(cos(y), 2));
    return z*pars[SIGMA] + pars[MEAN];
}

qreal MoyalDistr::func(qreal z)
{
    return exp( -0.5*(z + exp(-z)) )/sqrt(2*M_PI);
}

qreal MoyalDistr::getMean(const array1 &in)
{
    qreal out = 0;
    for (int i=0; i<in.size(); i++) {
        out += in.at(i)/in.size();
    }
    return out;
}

qreal MoyalDistr::getSigma(const array1 &in, qreal mean)
{
    qreal sum =0;
    if (mean == 0)
        mean = getMean(in);
    for (int i=0; i<in.size(); i++) {
        sum += pow(in.at(i) - mean, 2);
    }
    return sqrt(sum/(in.size()));
}

bool MoyalDistr::verify(const array1 &tmpPars)
{
//    return true;
    return tmpPars.at(SIGMA) > 0;
}
