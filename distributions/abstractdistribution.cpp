#include "abstractdistribution.h"

//int AbstractDistribution::count = 0;
RNG AbstractDistribution::rng;
AbstractDistribution::Estimation AbstractDistribution::est = AbstractDistribution::MD;

AbstractDistribution::AbstractDistribution()
{
    preDefined = false;
}

AbstractDistribution::AbstractDistribution(const array1 &in, const array1 *pars0)
{
    setEmpData(in);
    if (pars0 != 0) {
        preDefined = true;
        pars = *pars0;
    }
    else
        preDefined = false;
}

AbstractDistribution::~AbstractDistribution()
{

}

qreal AbstractDistribution::cdf(qreal x, const array1 &tmpPars)
{
    if (verify(tmpPars))
        return cdfVerified(x, tmpPars);

    return -1;
}

array1 AbstractDistribution::estimate(array1 in, array1 pars0)
{
    if (est == AbstractDistribution::MD)
        return findMDPars(in, pars0);
    else if (est == AbstractDistribution::MLH)
        return findMLHPars(in, pars0);
    else
        return pars0;
}

array1 AbstractDistribution::generateSeries(int n)
{
    array1 out;
    for (int i=0; i<n; i++) {
//        qDebug () << i;
        qreal tmp = rnd();
//        qDebug () << tmp;
        out << tmp/*rnd()*/;
    }

    return out;
}

array1 AbstractDistribution::generateSeries()
{
    return generateSeries(data.size());
}

array1 AbstractDistribution::findMDPars(array1 in, array1 pars0)
{
//    qDebug () << "MD-estimation" ;
    int N = 1000, M = 10;
//    M *= pars.size();
    qreal step = 4, beta=0.8, e=1E-3, basestep, stepVariance = 1;

    array1 x0, x1, ksi;
    array1 steps = pars0;
//    array2 matrix;
    x0 = pars0;
    x1 = x0;
    basestep = step;
    qreal min = fitCriteria(in, x0, omega_2, true);
    int i=0, j=0;
//    do {
//        array1 tmpArr (pars0.size());
//    }
//    while ();

    for (i=0; i<steps.size(); i++) {
        steps[i] *=.95;
    }

    for (i=0; i<N; i++) {

        bool chgCenter = false;

        for (j=0; j<M; j++) {
            ksi = vectorMult(randomUnitVector(x0/*, pars0*/), step);
            array1 tmp = vectorSum(x0, ksi);
//            array1 tmp = x0;
//            for (int k=0; k<x0.size(); k++) {
//                tmp[k] = ;
//            }
            if (!verify(tmp))
                continue;
            qreal S = fitCriteria(in, tmp, omega_2, true);
            if (S < min) {
                chgCenter = true;
                min = S;
                x1 = tmp;
                break;
            }
        }

        if (chgCenter) {
            x0 = x1;
            if (stepVariance < 2)
                stepVariance += 0.00 /*0.05*/;
            step = basestep*stepVariance;
        }
        else {
            step *= beta;
            basestep = step;
            stepVariance = 1;
        }

        if (step < e) {
            minOrMax = min;
            qDebug () << numGenerations();
            return x0;
        }

    }
    return x0;
}

array1 AbstractDistribution::findMLHPars(array1 in, array1 pars0)
{
//    qDebug () << "MLH-estimation" ;
    int N = 1000, M = 20;
    qreal step = 1, beta=0.95, e=1E-3, basestep, stepVariance = 1;
    array1 x0, x1, ksi;
    x0 = pars0;
    x1 = x0;
    basestep = step;

    for (int i=0; i<N; i++) {
        qreal max = likelihood(in, x0);
        bool chgCenter = false;

        for (int j=0; j<M; j++) {
            ksi = vectorMult(randomUnitVector(x0/*, pars0*/), step);
            array1 tmp = vectorSum(x0, ksi);
            if (!verify(tmp))
                continue;
            qreal S = likelihood(in, tmp);
            if (S > max) {
                chgCenter = true;
                max = S;
                x1 = tmp;
                break;
            }
        }

        if (chgCenter) {
            x0 = x1;
            if (stepVariance < 2)
                stepVariance += /*0.00*/ 0.05;
            step = basestep*stepVariance;
        }
        else {
            step *= beta;
            basestep = step;
            stepVariance = 1;
        }

        if (step < e) {

            minOrMax = max;
            return x0;
        }
//        qDebug () << i << step-e;
    }
//    pars0 = HookJeevs(in, x0);
    qDebug () << "> e!";
    return x0;
}

qreal AbstractDistribution::findPValue()
{
    qreal crit;
    return findPValue(data, crit);
}

qreal AbstractDistribution::findPValue(array1 in, qreal &criteria, int iters, Method method)
{

    if (method != Kolmogorov || method != Omega_2)
        method = omega_2;

    qreal pvalue=0;
    criteria = fitCriteria(in, method);


    for (int i=0; i<iters; i++) {
        array1 s = this->generateSeries();

//        pars << p;
        /*

        //! Двойной бутстреп
        array1 p = distr->getPars(s);

        //! меняем вручную параметры распределения
        distr->setParameter(p);
        //! Генерируем ещё одну выборку (вторая петля бутстрепа)
        s = distr->generateSeries();
        //! меняем параметры на первоначальные
        distr->setParameter(startPars);*/

        if (fitCriteria(s, method) > criteria)
            pvalue ++;
    }

    return pvalue*1.0/iters;

}

qreal AbstractDistribution::fitCriteria(array1 in, const array1 &tmpPars, Method method, bool chkPars)
{
    if (in.isEmpty())
        return 0;
    else {
        qSort(in);
    }
    int n = in.size();

//    if (name()=="DoubleNormalDistr" && chkPars) {
////        bool norm = (i==0 || i==2);
//        if (tmpPars[4] < 0 || tmpPars[4] > 1)
//            return 1E10;
//    }

//    for (int i=0; i<tmpPars.size(); i++) {
//        bool norm = (name()==QString("NormalDistr") && i==0) ||
//                    (name()==QString("DoubleNormalDistr") && (i==0 || i==2));
    //        if (tmpPars.at(i) <=0 && chkPars && !norm)
    //            return 1E10;
    //    }


    qreal cdfI;


    if (method == Kolmogorov) {
        qreal Dn=0, Dnp=0, Dnm=0;
        for (int i=1; i<=n; i++)
        {
            qreal d1, d2;
            cdfI = this->cdf(in.at(i-1), tmpPars);
            if (cdfI < 0)
                return 1E10;

            d1 = i*1.0/n - cdfI;
            d2 = cdfI - (i-1.0)/n;
            if (i == 0) {
                Dnp = d1;
                Dnm = d2;
            }

            if (d1 > Dnp)
                Dnp = d1;

            if (d2 > Dnm)
                Dnm = d2;


            if (Dnp > Dnm)
                Dn = Dnp;
            else
                Dn = Dnm;
        }

        return (6.0*n*Dn + 1.0) / (6.0*sqrt(n));
    }

    else if (method == Omega_2) {
        qreal summ = 0;


        for (int i=1; i<=n; i++) {
            cdfI = this->cdf(in.at(i-1), tmpPars);
            if (cdfI < 0)
                return 1E10;
            summ += (2.0*i - 1.0)/(2.0*n)*log(cdfI);
            summ += (1.0 - (2.0*i - 1.0)/(2.0*n))*log(1.0 - cdfI);
        }

        return -n - 2.0*summ;
    }

    else {
        qreal So = 1.0/(12*n);

        for (int i=1; i<=n; i++) {
            cdfI = this->cdf(in.at(i-1), tmpPars);
            if (cdfI < 0)
                return 1E10;
            qreal delta;
            delta = cdfI - (2.0*i - 1.0)/(2*n);
            So += pow(delta, 2);
        }


        return So;
    }

}

qreal AbstractDistribution::likelihood(const array1 &in, const array1 &tmpPars)
{
    qreal multi = 1.0;
//    multi = 0.0;
    for (int i=0; i<in.size(); i++)
    {
//        multi += log(pdf(in.at(i), tmpPars));
        multi *= pdf(in.at(i), tmpPars);
    }
    return multi;
}

void AbstractDistribution::setEmpData(const array1 &in)
{
    if (in.isEmpty())
        qDebug () << "Input array is empty!";

    data = in;
    sortedData = data;
    qSort(sortedData);
}

array1 AbstractDistribution::getPars(const array1 &input)
{
    array1 output;
    defineParameters(input, output);
    return output;
}

bool AbstractDistribution::setParameter(int id, qreal value)
{
    if (id < 0 || id >= pars.size())
        return false;

    pars[id] = value;
    return true;
}

bool AbstractDistribution::setParameter(const array1 &in)
{
    int size;
    if (in.size() <= pars.size())
        size = in.size();
    else
        size = pars.size();

    for (int i=0; i<size; i++) {
        pars[i] = in.at(i);
    }
    return true;
}

qreal AbstractDistribution::fitCriteria(array1 in, const array1 &tmpPars, qreal x, int ind)
{
    array1 tmp = tmpPars;
    tmp[ind] = x;
    return fitCriteria(in, tmp, omega_2, true);
}

void AbstractDistribution::find2DMin(array1 in, array1 &pars0, int ind)
{
    array1 x(4), y(4);
    const qreal goldR = 0.61803;
    const qreal eps = 1E-3;
    const qreal goldL = 1.0 - goldR;
    x[0] = pars0.at(ind)*0.1;
    x[3] = pars0.at(ind)*1.9;
    x[1] = x[0]+(x[3]-x[0])*goldL;
    x[2] = x[0]+(x[3]-x[0])*goldR;

    for (int i=0; i<4; i++) {
        y[i] = fitCriteria(in, pars0, x.at(i), ind);
    }

    while (fabs(x[3]-x[0]) > eps) {
        if (y[1] < y[2]) {
            //отбрасываем крайнюю правую точку
            x[3] = x[2];
            y[3] = y[2];
            x[2] = x[1];
            y[2] = y[1];
            x[1] = x[0]+(x[3]-x[0])*goldL;
            y[1] = fitCriteria(in, pars0, x.at(1), ind);
        }
        else if (y[2] < y[1]) {
            //отбрасываем крайнюю левую точку
            x[0] = x[1];
            y[0] = y[1];
            x[1] = x[2];
            y[1] = y[2];
            x[2] = x[0]+(x[3]-x[0])*goldR;
            y[2] = fitCriteria(in, pars0, x.at(2), ind);
        }
        else {
            qDebug () << x << y;
            return;
        }
    }
    pars0[ind] = 0.5*(x[3]+x[0]);
}

array1 AbstractDistribution::vectorMult(array1 in, qreal multiplicator)
{
    array1 out;
    for (int i=0; i<in.size(); i++) {
        out << in.at(i)*multiplicator;
    }
    return out;
}

array1 AbstractDistribution::vectorSum(array1 in1, array1 in2, qreal m1, qreal m2)
{
    array1 out;
    if (in1.size() != in2.size())
        return out;
    for (int i=0; i<in1.size(); i++) {
        out << m1*in1.at(i) + m2*in2.at(i);
    }
    return out;
}

array1 AbstractDistribution::randomUnitVector(array1 p)
{
    array1 out (p.size());
    qreal length = 0;
    for (int i=0; i<p.size(); i++) {
        qreal border = 1;
//        if (p.size() == 0)
//            border = 1;
//        else
        border = p.at(i)/2.0;
        out[i] = rng.uniform(-border, border);
        length += pow(out[i], 2);
    }
    length = sqrt (length);
    for (int i=0; i<p.size(); i++) {
        out[i] /= length;
    }
    return out;
}

qreal AbstractDistribution::explore(array1 in, array1 h, array1 x1, array1 &x2)
{
    int n = x1.size();
    qreal y1, y2;
    x2 = x1;

    y1 = fitCriteria(in, x1);
    for (int i=0; i<n; i++) {
        y2 = fitCriteria(in, x2, x2[i]+h[i], i);
        if (y2 < y1) {
            x2[i] = x1[i] + h[i];
            continue;
        }
        y2 = fitCriteria(in, x2, x2[i]-h[i], i);
        if (y2 < y1) {
            x2[i] = x2[i] - h[i];
            continue;
        }
    }
    return fitCriteria(in, x2);
}

array1 AbstractDistribution::HookJeevs(array1 in, array1 pars0)
{
//    qDebug () << pars0;
//    pars0[0] = 1.80569;
//    pars0[1] = 0.580693;
    int n = pars0.size();
    array1 h (n);
    qreal y1, y2, e=1E-3, chk=0;
    array1 x1 (n), x2 (n), x3 (n), P (n);
    x1 = pars0;
//    x2 = x1;
    for (int i=0; i<n; i++) {
        h[i] = pars0[i]*0.95;
        h[i] = .5;
        chk += pow (h[i], 2);
    }

    chk = sqrt (chk);
    int N = 20;
    int NN = 0;

    while (chk > e) {
        //Исследование
        /*y1 = fitCriteria(in, x1);
        for (int i=0; i<n; i++) {
            y2 = fitCriteria(in, x1, x1[i]+h, i);
            if (y2 < y1) {
                x2[i] += h;
                ptrn[i] = 1;
                continue;
            }
            y2 = fitCriteria(in, x1, x1[i]-h, i);
            if (y2 < y1) {
                x2[i] -= h;
                ptrn[i] = -1;
                continue;
            }
        }*/


        y1 = explore(in, h, x1, x2);

        if (x1 == x2) {
            h = vectorMult(h, .1);
            chk = 0;
            for (int i=0; i<n; i++) {
                chk += pow(h[i], 2);
            }
            chk = sqrt (chk);
            continue;
        }

        //Поиск по образцу
        do {
            for (int i=0; i<n; i++) {
                P[i] = x1[i] + 2*(x2[i] - x1[i]);
            }

            y2 = explore (in, h, P, x3);
            if (y2 < y1) {
                x2 = x3;
            }
            else {
                x1 = x2;
                break;
            }
        } while (true);

        chk = 0;
        for (int i=0; i<n; i++) {
            chk += pow(h[i], 2);
        }
        chk = sqrt (chk);
        NN++;
        if (NN >= N)
            break;
    }
    return x1;
}
