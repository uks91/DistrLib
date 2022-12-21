#ifndef MOYALDISTR_H
#define MOYALDISTR_H

#include "abstractdistribution.h"

class MoyalDistr : public AbstractDistribution
{
public:

    enum {MEAN, SIGMA};

    MoyalDistr();
    MoyalDistr(const array1 &in, const array1 *pars0=0);

    qreal cdfVerified (qreal x, const array1 &tmpPars);
//    qreal optimum (qreal mean, qreal sigma, array1 in);
    void defineParameters(const array1 &d, array1 &pars);
    QString name () {return "MoyalDistr";}
    qreal pdf (qreal x, const array1 &tmpPars);
    qreal rnd();

private:
    qreal func (qreal z);
    qreal getMean (const array1 &in);
    qreal getSigma (const array1 &in, qreal mean=0);
    bool  verify(const array1 &tmpPars);
};

#endif // MOYALDISTR_H
