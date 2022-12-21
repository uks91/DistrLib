#ifndef NORMALDISTR_H
#define NORMALDISTR_H

#include "abstractdistribution.h"

class NormalDistr : public AbstractDistribution
{

public:
    enum {MEAN, SIGMA};

    NormalDistr(/*bool isMirrored=false, qreal shift=0*/);
    NormalDistr(const array1 &in, const array1 *pars0=0);

    qreal cdfVerified (qreal x, const array1 &tmpPars);
    qreal optimum (qreal mean, qreal sigma, array1 in);
    void defineParameters(const array1 &d, array1 &pars);
    QString name () {return "NormalDistr";}
    qreal pdf (qreal x, const array1 &tmpPars);
    qreal rnd();

protected:
    qreal getMean (const array1 &in);
    qreal getSigma (const array1 &in, qreal mean=0);
    qreal getMean ();
    qreal getSigma ();
    bool  verify(const array1 &tmpPars);

};


#endif // NORMALDISTR_H
