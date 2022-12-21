#ifndef LOGNORMALDISTR_H
#define LOGNORMALDISTR_H

#include "normaldistr.h"

class LognormalDistr : public NormalDistr
{
public:
    LognormalDistr();

//    LognormalDistr(/*bool isMirrored=false, qreal shift=0*/);
    LognormalDistr(const array1 &in, const array1 *pars0=0);

    qreal cdfVerified (qreal x, const array1 &tmpPars);
    QString name () {return "LognormalDistr";}
//    void defineParameters();
//    qreal optimum (qreal mean, qreal sigma, array1 in);
    void defineParameters(const array1 &d, array1 &pars);
    qreal pdf (qreal x, const array1 &tmpPars);
    qreal rnd();

private:
    array1 raw_data;

//    void setLogData ();
    array1 logData(const array1 &in);
    bool   verify(const array1 &tmpPars);

};

#endif // LOGNORMALDISTR_H
