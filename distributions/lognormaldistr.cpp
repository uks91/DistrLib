#include "lognormaldistr.h"

LognormalDistr::LognormalDistr (/*bool isMirrored, qreal shift*/) : NormalDistr()
{
    pars.resize(2);
}


LognormalDistr::LognormalDistr(const array1 &in, const array1 *pars0): NormalDistr()
{
    setEmpData(in);
    raw_data = in;
    if (pars0 != 0) {
        pars = *pars0;
        preDefined = true;
    }
    else
        preDefined = false;
//    data = logData(in);

    AbstractDistribution::defineParameters();
}

qreal LognormalDistr::cdfVerified(qreal x, const array1 &tmpPars)
{
    if (x<=0)
        return 0;
//    qreal tmp = NormalDistr::cdf(log(x), tmpPars);
    return NormalDistr::cdfVerified(log(x), tmpPars);
}

//void LognormalDistr::defineParameters()
//{

//}

/*qreal LognormalDistr::optimum (qreal mean, qreal sigma, array1 in)
{
    LognormalDistr *distr = new LognormalDistr();
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
}*/

void LognormalDistr::defineParameters(const array1 &d, array1 &pars)
{
    pars.resize(2);
    array1 dd = logData(d);
    pars[MEAN]  = getMean(dd);
    pars[SIGMA] = getSigma(dd);

    pars = estimate(d, pars);
    //    pars=findMDPars(d, pars);

    return;

    qreal min = 1E9;
    qreal mean, sigma;
    qreal mean0, sigma0;
    mean0 = pars[MEAN];
    sigma0 = pars[SIGMA];
    for (mean=mean0*0.5; mean <=1.5*mean0; mean += 0.01*mean0) {
        for (sigma=0.5*sigma0; sigma <=1.5*sigma0; sigma += 0.01*sigma0) {
            qreal tmp = optimum(mean, sigma, d);
            if (tmp < min) {
                pars[MEAN] = mean;
                pars[SIGMA] = sigma;
                min = tmp;
            }
        }
    }
}

qreal LognormalDistr::pdf(qreal x, const array1 &tmpPars)
{
    if (x<0)
        return 0;
    return NormalDistr::pdf(log(x), tmpPars)/x;
}

qreal LognormalDistr::rnd()
{
//    return (NormalDistr::rnd());
    return exp(NormalDistr::rnd());
}

//void LognormalDistr::setLogData()
//{
////    logs = data;
////    for (int i=0; i<logs.size(); )
//}

array1 LognormalDistr::logData(const array1 &in)
{
    array1 out = in;
    for (int i=0; i<out.size(); i++) {
        out[i] = log(in.at(i));
//        out[i] = (in.at(i));
    }

    return out;
}

bool LognormalDistr::verify(const array1 &tmpPars)
{
    return tmpPars.at(SIGMA) > 0;
}
