#ifndef DOUBLENORMALDISTR_H
#define DOUBLENORMALDISTR_H

#include "abstractdistribution.h"
#include "normaldistr.h"

class DoubleNormalDistr : public AbstractDistribution
{
public:    
    enum {MEAN1, SIGMA1, MEAN2, SIGMA2, A};

    DoubleNormalDistr(/*bool isMirrored=false, qreal shift=0*/);
    DoubleNormalDistr(const array1 &in, const array1 *pars0=0);

    qreal cdfVerified (qreal x, const array1 &tmpPars);
    qreal optimum (qreal mean, qreal sigma, array1 in);
    void defineParameters(const array1 &d, array1 &pars);
    QString name () {return "DoubleNormalDistr";}
    qreal pdf (qreal x, const array1 &tmpPars);
    qreal rnd();

protected:
    qreal getMean (const array1 &in);
    qreal getSigma (const array1 &in, qreal mean=0);
    qreal getMean ();
    qreal getSigma ();
    NormalDistr n1;
    NormalDistr n2;
    bool prepareNormals(const array1 &in);
    bool verify(const array1 &tmpPars);

};

#endif // DOUBLENORMALDISTR_H
