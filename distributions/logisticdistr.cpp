#include "logisticdistr.h"

LogisticDistr::LogisticDistr() : AbstractDistribution()
{
    pars.resize(2);

    //    AbstractDistribution::defineParameters();
}

LogisticDistr::LogisticDistr(const array1 &in, const array1 *pars0) : AbstractDistribution(in, pars0)
{
    pars.resize(2);

    AbstractDistribution::defineParameters();
}

qreal LogisticDistr::cdfVerified(qreal x, const array1 &tmpPars)
{
    qreal z = (x - tmpPars[ALPHA])/tmpPars[BETA];
    return 1.0 / (1.0 + exp(-z));
}

void LogisticDistr::defineParameters(const array1 &d, array1 &pars)
{
    if (pars.size() !=2)
        pars.resize(2);
    pars[ALPHA] = getMean(d);
    pars[BETA] = getSigma(d, pars[ALPHA]);

//    qDebug () << "CDF: " << cdf (0.5, pars);

    pars = estimate(d, pars);

    //    pars = findMDPars(d, pars);
}

qreal LogisticDistr::pdf(qreal x, const array1 &tmpPars)
{
    qreal z = (x - tmpPars[ALPHA])/tmpPars[BETA];
    return exp(z)*pow(1.0 + exp(z), -2)/tmpPars[BETA];
}

qreal LogisticDistr::rnd()
{
    qreal y = rng.uniform(0, 1);
    return pars[ALPHA] + pars[BETA]*log(y/(1.0 - y));
}

qreal LogisticDistr::getMean(const array1 &in)
{
    qreal out = 0;
    for (int i=0; i<in.size(); i++) {
        out += in.at(i)/in.size();
    }
    return out;
}

qreal LogisticDistr::getSigma(const array1 &in, qreal mean)
{
    qreal sum =0;
    if (mean == 0)
        mean = getMean(in);
    for (int i=0; i<in.size(); i++) {
        sum += pow(in.at(i) - mean, 2);
    }
    return sqrt(sum/(in.size()));
}

bool LogisticDistr::verify(const array1 &tmpPars)
{
    return tmpPars.at(BETA) > 0;
}
