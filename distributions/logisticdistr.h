#ifndef LOGISTICDISTR_H
#define LOGISTICDISTR_H

#include "abstractdistribution.h"

class LogisticDistr : public AbstractDistribution
{
public:

    enum {ALPHA, BETA};

    LogisticDistr();
    LogisticDistr(const array1 &in, const array1 *pars0=0);

    qreal cdfVerified (qreal x, const array1 &tmpPars);
    void defineParameters(const array1 &d, array1 &pars);
    QString name () {return "LogisticDistr";}
    qreal pdf (qreal x, const array1 &tmpPars);
    qreal rnd();

private:
    qreal getMean (const array1 &in);
    qreal getSigma (const array1 &in, qreal mean=0);
    bool  verify(const array1 &tmpPars);
};

#endif // LOGISTICDISTR_H
