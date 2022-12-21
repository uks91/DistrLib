#include "doublenormaldistr.h"


DoubleNormalDistr::DoubleNormalDistr (/*bool isMirrored, qreal shift*/) : AbstractDistribution()
{
    pars.resize(5);
}

DoubleNormalDistr::DoubleNormalDistr(const array1 &in, const array1 *pars0): AbstractDistribution(in)
{
    pars.resize(5);

    AbstractDistribution::defineParameters();
}

qreal DoubleNormalDistr::cdfVerified(qreal x, const array1 &tmpPars)
{
    prepareNormals(tmpPars);

    return tmpPars[A]*n1.AbstractDistribution::cdf(x) + (1.0-tmpPars[A])*n2.AbstractDistribution::cdf(x);

}

void DoubleNormalDistr::defineParameters(const array1 &d, array1 &pars)
{
    pars.resize(5);
    n1 = NormalDistr(d/*, pars*/);
    n2 = n1;
    qreal mean   = n1.getPars().at(NormalDistr::MEAN);
    qreal sigma  = n1.getPars().at(NormalDistr::SIGMA);
    qreal m = 0.5;
    pars[MEAN1]  = mean*rng.uniform(-m, m);
    pars[SIGMA1] = sigma*rng.uniform(-m, m);
    pars[MEAN2]  = mean*rng.uniform(-m, m);
    pars[SIGMA2] = sigma*rng.uniform(-m, m);
    pars[A]      = rng.uniform(0, 1);

    //2.84041, 1.52895, 1.67582, 0.513171, 0.185635 //0.044 (0.913)

    qreal min0=1E10, max0=0;
    array1 tmpPars = pars;
//    min = 1E9;
    int numIters = 20;
//    numIters = 1;

    for (int i=0; i<numIters; i++) {
        qreal mean   = n1.getPars().at(NormalDistr::MEAN);
        qreal sigma  = n1.getPars().at(NormalDistr::SIGMA);
        qreal m = 0.5;
        tmpPars[MEAN1]  = mean*rng.uniform(-m, m);
        tmpPars[SIGMA1] = sigma*rng.uniform(-m, m);
        tmpPars[MEAN2]  = mean*rng.uniform(-m, m);
        tmpPars[SIGMA2] = sigma*rng.uniform(-m, m);
        tmpPars[A]      = rng.uniform(0, 1);

        tmpPars  = estimate(d, tmpPars);
//        qDebug () << "Max: " << minOrMax ;
        if ((minOrMax < min0) && est==AbstractDistribution::MD) {
            pars = tmpPars;
            min0 = minOrMax;
//            qDebug () << pars << min;
        }

        if ((minOrMax > max0) && est==AbstractDistribution::MLH /*&& minOrMax != -INFINITY*/) {
            pars = tmpPars;
            max0 = minOrMax;
            qDebug () << pars << minOrMax;
        }
    }

//    pars = findMDPars(d, pars);
//    pars = findMDPars(d, pars);
//    pars = findMDPars(d, pars);

//    qDebug () << pars;
}

qreal DoubleNormalDistr::getMean()
{
    return getMean(data);
}

qreal DoubleNormalDistr::getMean(const array1 &in)
{
    qreal out = 0;
    for (int i=0; i<in.size(); i++) {
        out += in.at(i)/in.size();
    }
    return out;
}

qreal DoubleNormalDistr::getSigma()
{
    return getSigma(data);
}

qreal DoubleNormalDistr::getSigma(const array1 &in, qreal mean)
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

qreal DoubleNormalDistr::pdf(qreal x, const array1 &tmpPars)
{
    if (prepareNormals(tmpPars))
        return tmpPars[A]*n1.AbstractDistribution::pdf(x) + (1.0-tmpPars[A])*n2.AbstractDistribution::pdf(x);
    else
        return -1;
}

qreal DoubleNormalDistr::rnd()
{
    prepareNormals(pars);
//    return pars[A]*n1.rnd() + (1.0-pars[A])*n2.rnd();

    const qreal e=1e-4;

    qreal res = rng.uniform(0, 1);

    array1 x (3);
    array1 u (2);

    x[1] = pars[MEAN1];  //Начальное приближение
    x[0] = x[1]*0.85;  //Второе начальное приближение

    do {

        u[0] = cdf(x[0]) - res;
        u[1] = cdf(x[1]) - res;

        x[2] = x[1] - u[1]*(x[1]-x[0])/(u[1]-u[0]);
        x[0] = x[1];
        x[1] = x[2];
    }
    while (fabs(/*res - */u[1])>e);

    return x[1];
}

bool DoubleNormalDistr::prepareNormals(const array1 &in)
{
    if (!verify(in))
        return false;

    n1.setParameter(NormalDistr::MEAN,  in.at(MEAN1));
    n1.setParameter(NormalDistr::SIGMA, in.at(SIGMA1));
    n2.setParameter(NormalDistr::MEAN,  in.at(MEAN2));
    n2.setParameter(NormalDistr::SIGMA, in.at(SIGMA2));

    return true;
}

bool DoubleNormalDistr::verify(const array1 &tmpPars)
{
    return !(tmpPars.at(SIGMA1) <= 0 ||
            tmpPars.at(SIGMA2) <= 0 ||
            tmpPars.at(A) <=0 ||
            tmpPars.at(A) >= 1);
}
