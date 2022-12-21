#ifndef WEIBULLDISTR_H
#define WEIBULLDISTR_H

#include "abstractdistribution.h"

class WeibullDistr : public virtual AbstractDistribution
{


public:
    enum {SHAPE, SCALE, SHIFT};
    WeibullDistr(bool isMirrored=false, qreal shift=0);
    WeibullDistr(const array1 &in, const array1 *pars0=0);

    qreal optimum (qreal scale, qreal shape, array1 in);
    qreal cdfVerified (qreal x, const array1 &tmpPars);
//    void defineShapeFromScale ();
    qreal defineScaleFromShape (qreal shape, const array1 &in);
    void defineScaleFromShape () {pars[SCALE] = defineScaleFromShape(pars[SHAPE], data);}
    void defineParameters(const array1 &in, array1 &pars);
    QString name () {return "WeibullDistr";}
    qreal pdf (qreal x, const array1 &tmpPars);
    qreal rnd();
//    void setEmpData(const array1 &in);

private:
    bool isMirrored;
    void getMirroredData ();
    bool verify(const array1 &tmpPars);

    array1 mirroredData;

};

#endif // WEIBULLDISTR_H
