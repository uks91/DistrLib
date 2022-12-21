#include "extremevaluedistr.h"


ExtremeValueDistr::ExtremeValueDistr() : AbstractDistribution()
{
    pars.resize(2);

//    AbstractDistribution::defineParameters();
}

ExtremeValueDistr::ExtremeValueDistr(const array1 &in, const array1 *pars0) : AbstractDistribution(in, pars0)
{
    pars.resize(2);

    AbstractDistribution::defineParameters();
}

qreal ExtremeValueDistr::cdfVerified(qreal x, const array1 &tmpPars)
{
    qreal z = (x - tmpPars[ALPHA])/tmpPars[BETA];
    return exp(-exp(-z));
}

void ExtremeValueDistr::defineParameters(const array1 &d, array1 &pars)
{
    if (pars.size() !=2)
        pars.resize(2);
    pars[ALPHA] = getMean(d);
    pars[BETA] = getSigma(d, pars[ALPHA]);

    pars = estimate(d, pars);
//    pars = findMDPars(d, pars);
}

qreal ExtremeValueDistr::pdf(qreal x, const array1 &tmpPars)
{
    qreal z = (x - tmpPars[ALPHA])/tmpPars[BETA];
    return exp(-z-exp(-z))/tmpPars[BETA];
}

qreal ExtremeValueDistr::rnd()
{
    qreal y = rng.uniform(0, 1);
    return pars[ALPHA] - pars[BETA]*log(-log(y));
}

qreal ExtremeValueDistr::getMean(const array1 &in)
{
    qreal out = 0;
    for (int i=0; i<in.size(); i++) {
        out += in.at(i)/in.size();
    }
    return out;
}

qreal ExtremeValueDistr::getSigma(const array1 &in, qreal mean)
{
    qreal sum =0;
    if (mean == 0)
        mean = getMean(in);
    for (int i=0; i<in.size(); i++) {
        sum += pow(in.at(i) - mean, 2);
    }
    return sqrt(sum/(in.size()));
}

bool ExtremeValueDistr::verify(const array1 &tmpPars)
{
    return !(tmpPars.at(BETA) <= 0);
}
