#include "empiricaldistr.h"

EmpiricalDistr::EmpiricalDistr() : AbstractDistribution()
{
}

EmpiricalDistr::EmpiricalDistr(const array1 &in) : AbstractDistribution(in)
{
    qSort (data);
    AbstractDistribution::defineParameters();
}

void EmpiricalDistr::defineParameters(const array1 &data, array1 &/*pars*/)
{
    if (data.size() < 2) {
        qDebug () << "Too small" << data;
        return;
    }

    array1 in = data;

    qSort (in);
    expCPD.clear();
    array1 x, p;

    x << in.at(0);
    p << 0;

    //qDebug () << x.last() << p.last();

    for (int i=1; i<in.size(); i++) {
        if (in.at(i) != x.last()) {
            x << in.at(i);
            p << i*1.0/data.size();

            //qDebug () << x.at(x.size()-2)+1E-4 << p.last();
            //qDebug () << x.last() << p.last();

        }
    }

    qreal x0, x1, p0, p1;

    x0 = 1.5*x.at(0) - .5*x.at(1);
    p0 = 0; //return;

    array1 tmp;
    qreal k, b;


    for (int i=0; i<x.size()-1; i++) {

        x1 = (x.at(i+1) + x.at(i)) / 2.0;
        p1 = p.at(i+1);
        k = (p1 - p0)/(x1 - x0);
        b = p1 - k*x1;

        tmp << x0 << x1 << p0 << p1 << k << b;
        //qDebug () << x0 << x1 << p0 << p1 << k << b;

        x0 = x1;
        p0 = p1;
        expCPD << tmp;
        tmp.clear();
    }

    x1 = x.last() + (x.last() - x.at(x.size()-2))/2;
    tmp << x0 << x1 << p0 << 1 << k << b;
    expCPD << tmp;

//    qDebug () << expCPD;
}

qreal EmpiricalDistr::quantile (qreal p)
{
    qreal k, b;
    int i=0;
    for (i=0; i<expCPD.size(); i++) {
        if (expCPD.at(i).at(2) <= p && p <= expCPD.at(i).at(3))
            break;
    }

    k = expCPD.at(i).at(4);
    b = expCPD.at(i).at(5);

    return (p - b)/k;
}

qreal EmpiricalDistr::rnd()
{
    int p = rng.uniform(0, data.size()-1);
//    qDebug () << p << data.size();
    return data.at(p);
    return quantile(rng.uniform(0,1));
}

int EmpiricalDistr::realToInt(qreal x)
{
    qreal m = 0.5*(floor(x) + ceil(x));
    if (x >= m)
        return ceil(x);
    else
        return floor (x);
}
