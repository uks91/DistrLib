#include "normaldistr.h"

NormalDistr::NormalDistr (/*bool isMirrored, qreal shift*/) : AbstractDistribution()
{
    pars.resize(2);
}

NormalDistr::NormalDistr(const array1 &in, const array1 *pars0): AbstractDistribution(in)
{
    pars.resize(2);

    AbstractDistribution::defineParameters();
}

qreal NormalDistr::cdfVerified(qreal x, const array1 &tmpPars)
{
//    qreal tmp;
//    tmp = x - tmpPars[MEAN];
//    tmp = tmpPars[SIGMA]*sqrt(2);
//    tmp = (x - tmpPars[MEAN]) / (tmpPars[SIGMA]*sqrt(2));
//    tmp = erf (tmp);
//    tmp += 1;
//    tmp *= 0.5;
    return 0.5*(1 + erf( (x - tmpPars[MEAN]) / (tmpPars[SIGMA]*sqrt(2)) ));
}

qreal NormalDistr::optimum (qreal mean, qreal sigma, array1 in)
{
    NormalDistr *distr = new NormalDistr();
    int n = in.size();
    qreal So = 1.0/(12*n);


    distr->setParameter(MEAN, mean);
    distr->setParameter(SIGMA, sigma);

    qSort(in);

    for (int i=1; i<=n; i++) {
        qreal delta;
        delta = distr->AbstractDistribution::cdf(in.at(i-1)) - (2.0*i - 1.0)/(2*n);
        So += pow(delta, 2);
    }

    delete distr;

    return So;
}

void NormalDistr::defineParameters(const array1 &d, array1 &pars)
{
    pars.resize(2);
    pars[MEAN]  = getMean(d);
    pars[SIGMA] = getSigma(d, pars[MEAN]);

    pars = estimate(d, pars); return;

    pars = findMDPars(d, pars);

    return;

    qreal min = 1E9;
    qreal mean, sigma;
    qreal mean0, sigma0;
    mean0 = pars[MEAN];
    sigma0 = pars[SIGMA];
    for (mean=mean0*0.5; mean <=1.5*mean0; mean += 0.01*mean0) {
//        sigma = getSigma(d,mean);
        for (sigma=0.5*sigma0; sigma <=1.5*sigma0; sigma += 0.01*sigma0) {
            qreal tmp = optimum(mean, sigma, d);
            if (tmp < min) {
                pars[MEAN] = mean;
                pars[SIGMA] = sigma;
                min = tmp;
            }
        }
    }

//    qDebug () << pars;
}

qreal NormalDistr::getMean()
{
    return getMean(data);
}

qreal NormalDistr::getMean(const array1 &in)
{
    qreal out = 0;
    for (int i=0; i<in.size(); i++) {
        out += in.at(i)/in.size();
    }
    return out;
}

qreal NormalDistr::getSigma()
{
    return getSigma(data);
}

qreal NormalDistr::getSigma(const array1 &in, qreal mean)
{
    qreal sum =0;
    if (mean == 0)
        mean = getMean(in);
    for (int i=0; i<in.size(); i++) {
        sum += pow(in.at(i) - mean, 2);
    }
    return sqrt(sum/(in.size()));
}


//qreal NormalDistr::rnd()
//{
//    qreal y = rng.uniform(0,1);

//}

qreal NormalDistr::pdf(qreal x, const array1 &tmpPars)
{
    qreal m = tmpPars[MEAN];
    qreal s = tmpPars[SIGMA];
    return exp(-.5*pow((x-m)/s, 2))/(s*sqrt(2*M_PI));
}

qreal NormalDistr::rnd()
{
    if (pars[SIGMA] == 0)
        return pars[MEAN];

    double x, y, r2;

    do
    {
        x = rng.uniform(-1, 1);
        y = rng.uniform(-1, 1);
        r2 = x*x + y*y;
    }
    while (r2 > 1.0 || r2 == 0);

    return pars[MEAN] + pars[SIGMA] * y * sqrt (-2.0 * log (r2) / r2);
}

bool NormalDistr::verify(const array1 &tmpPars)
{
    return tmpPars.at(SIGMA) > 0;
}
