#include "weibulldistr.h"

WeibullDistr::WeibullDistr (bool isMirrored, qreal shift) : AbstractDistribution()
{
    pars.resize(2);
    this->isMirrored = isMirrored;
//    pars[SHIFT] = shift;
}

WeibullDistr::WeibullDistr(const array1 &in, const array1 *pars0): AbstractDistribution(in, pars0)
{
    pars.resize(2);

//    if (in.isEmpty())
//        qDebug () << "Array is empty!";

//    pars[SHIFT] = shift;


//    getMirroredData();
    AbstractDistribution::defineParameters();
}

qreal WeibullEq(qreal shape, array1 in)
{
    //! Используется в системе уравнений при определении параметров распределения Вейбулла
    //! in - выборка
    //! shape - аргумент функции (который ищется при помощи getRoot)


    qreal sum1=0, sum2=0, sum3 = 0;
    int n = in.size();

    for (int i=0; i<n; i++) {
        sum1 += pow(in.at(i), shape)*log(in.at(i));
        sum2 += pow(in.at(i), shape);
        sum3 += log(in.at(i));
    }

    return sum1/sum2 - 1/shape - sum3/n;
}

qreal WeibullDistr::optimum (qreal scale, qreal shape, array1 in)
{
//    array1 in = data;
    WeibullDistr *distr = new WeibullDistr();
//    qreal mean = 0;
    int n = in.size();
    qreal So = 1.0/(12*n);

//    for (int i=0; i<n; i++)
//        mean += data.at(i)/n;

    distr->setParameter(SCALE, scale);
    distr->setParameter(SHAPE, shape);

//    if (in.isEmpty())
//        in = distr->getSortedData();
//    else {
//        qSort(in);
//    }
    qSort(in);

    for (int i=1; i<=n; i++) {
        qreal delta;
        delta = distr->AbstractDistribution::cdf(in.at(i-1)) - (2.0*i - 1.0)/(2*n);
        So += pow(delta, 2);
    }

//    qDebug () << So;

    delete distr;

    return So;
}

qreal WeibullDistr::cdfVerified(qreal x, const array1 &tmpPars)
{
//    qDebug () << "Weibull!";
//    if (!isMirrored) {
    array1 tmp = tmpPars;
    tmp << 0;
//    if (x <= tmp[SHIFT])
//        return 0;

//    tmp[SHIFT] = 0;

    return 1 - exp(-pow((x - tmp[SHIFT])/tmp[SCALE], tmp[SHAPE]));
//    }

//    if (isMirrored && x >= pars[SHIFT])
//        return 0;
    return 0;

}




void WeibullDistr::defineParameters(const array1 &in, array1 &pars)
{
//    array1 in;
//    if (isMirrored)
//        in = mirroredData;
//    else
//        in = data;

    pars.resize(2);

#define MD_WEIBULL

#ifdef MD_WEIBULL

//    qreal x1, x2;
//    x1 = .001;
//    qreal step=1;
//    qreal e=1E-4;
//    short int dir = 1;
//    qreal y1 = optimum(defineScaleFromShape(x1, in), x1, in);
//    qreal y2;
//    while (step>e) {
//        x2 = x1 + dir*step;
//        y2 = optimum(defineScaleFromShape(x2, in), x2, in);
//        if (y2 > y1) {
//            step *= .5;
//            dir *= -1;
//        }
//        x1 = x2;
//        y1 = y2;
//    }

//    pars[SCALE] = defineScaleFromShape(x2, in);
//    pars[SHAPE] = x2;

//    return;

    qreal shape0 = 2;
    shape0 = getRoot(WeibullEq, in, shape0); //метод максимального правдободобия

    qreal scale0 = 0;
    for (int i=0; i<in.size(); i++) {
        scale0 += pow(in.at(i), shape0)/in.size();
    }

    scale0 = pow(scale0, 1.0/shape0);

    pars[SCALE] = scale0;
    pars[SHAPE] = shape0;
    pars = estimate(in, pars); return;

    pars = findMDPars(in, pars); return;

    qreal shape = 1.8;
    qreal scale;
    qreal min = 1E9;

    qreal percent = 0.3;

    for (shape=(1-percent)*shape0; shape<=(1+percent)*shape0; shape += .01*shape0) {
//        scale = 0;
//        for (int i=0; i<in.size(); i++) {
//            scale += pow(in.at(i), shape)/in.size();
//        }

        //        scale = pow(scale, 1.0/shape);
        for (scale=(1-percent)*scale0; scale<=(1+percent)*scale0; scale += 0.01*scale0) {

            qreal p = optimum(scale, shape, in);
            if (min > p) {
                min = p;
                pars[SCALE] = scale;
                pars[SHAPE] = shape;
            }
        }
    }

    return;
#else

    pars[SHAPE] = 2;
    pars[SHAPE] = getRoot(WeibullEq, in, pars[SHAPE]); //метод максимального правдободобия

    pars[SCALE] = 0;
    for (int i=0; i<in.size(); i++) {
        pars[SCALE] += pow(in.at(i), pars[SHAPE])/in.size();
    }

    pars[SCALE] = pow(pars[SCALE], 1.0/pars[SHAPE]);
#endif
}

qreal WeibullDistr::defineScaleFromShape(qreal shape, const array1 &in)
{
    if (in.size() < 2)
        return 0;

    qreal scale;

    scale = 0;
    for (int i=0; i<in.size(); i++) {
        scale += pow(in.at(i), shape)/in.size();
    }

    scale = pow(scale, 1.0/shape);

    return scale;
}

void WeibullDistr::getMirroredData()
{
    if (!isMirrored)
        return;

//    mirroredData.clear();
//    for (int i=0; i<data.size(); i++) {
//        mirroredData << pars[SHIFT] - data.at(i);
//    }
}

qreal WeibullDistr::pdf(qreal x, const array1 &tmpPars)
{
    if (x<0)
        return 0;
    qreal b = tmpPars[SHAPE];
    qreal a = tmpPars[SCALE];
    return b/a*pow(x/a, b-1)*exp(-pow(x/a, b));
}

qreal WeibullDistr::rnd()
{
    qreal y = rng.uniform(0,1);
    qreal x = pars[SCALE]*pow(-log(1-y), 1.0/pars[SHAPE]);

    return x;
}

bool WeibullDistr::verify(const array1 &tmpPars)
{
    return tmpPars.at(SHAPE) > 0 && tmpPars.at(SCALE) > 0 ;
}
