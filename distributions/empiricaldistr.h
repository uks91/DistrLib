#ifndef EMPIRICALDISTR_H
#define EMPIRICALDISTR_H

#include "abstractdistribution.h"

class EmpiricalDistr : public AbstractDistribution
{
public:
    EmpiricalDistr();
    EmpiricalDistr(const array1 &in);

    qreal cdfVerified(qreal /*x*/, const array1 & /*pars*/) {return 0;}
    void defineParameters(const array1 &data, array1 &);
    QString name () {return "EmpiricalDistr";}
    qreal quantile (qreal p);
    qreal pdf (qreal /*x*/, const array1 &/*tmpPars*/) {return 0;}
    qreal rnd();

private:
    QVector<array1 > expCPD;
    int realToInt (qreal x);
    bool verify(const array1 &/*tmpPars*/) {return true;}

};

#endif // EMPIRICALDISTR_H
