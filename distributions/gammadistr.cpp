#include "gammadistr.h"
#include <specialfunctions.h>

//qreal psi (int k)
//{
//    qreal out = 0;
//    for (int i=0; i<200; i++) {
//        out += pow (i, -k);
//    }

//    return out;
//}

//qreal DLogGamma (qreal x)
//{
//    qreal out;
//    out = 1/(x+1) + 1-GammaDistr::Euler;

//    for (int k=2; k<100; k++) {
//        out += pow(-1, k)*(psi(k) - 1)*pow(x, k-1);
//    }
//    return out;
//}


qreal GammaEq(qreal beta, array1 in)
{
    //beta - scale
    //alpha - shape

    qreal alpha, sum = 0, logsum = 0;
    int n = in.size();
    for (int i=0; i<n; i++) {
        sum += in.at(i);
        logsum += log(in.at(i));
    }

    alpha  = sum/(n*beta)/* - 1*/;
    return logsum - n*log(beta) - n*alglib::psi(alpha);
}

GammaDistr::GammaDistr (qreal shift, bool isMirrored) : AbstractDistribution()
{
    pars.resize(2);

//    pars[ALPHA] = 1.3;
//    pars[BETA]  = 4;
//    pars[K]     = 2.3;

//    if (isMirrored)
    this->isMirrored = isMirrored;


//    pars[SHIFT] = shift;
}

GammaDistr::GammaDistr(const array1 &in, const array1 *pars0): AbstractDistribution(in, pars0)
{
    pars.resize(2);

//    if (isMirrored)
    this->isMirrored = isMirrored;

    //    if (in.isEmpty())
    //        qDebug () << "Array is empty!";

//    pars[SHIFT] = shift;

    getMirroredData();
    AbstractDistribution::defineParameters();
}

qreal GammaDistr::cdfVerified(qreal x, const array1 &tmpPars)
{
    return cpd(x, tmpPars);
    //    return gammaFunc(pars[ALPHA]+1, x/pars[BETA])/GammaFunc(pars[ALPHA]+1);

}

qreal /*GammaDistr::*/cpd(qreal x, array1 p)
{
    //    alpha--;
    if (x<0)
        return 0;
    qreal a,b;

    a = p[GammaDistr::ALPHA];
    b = p[GammaDistr::BETA];
    return GammaDistr::gammaFunc(a, x/b)/GammaDistr::GammaFunc(a);

}

qreal GammaDistr::optimum (qreal beta, array1 in, qreal alpha)
{
//    array1 in = data;
    GammaDistr *distr = new GammaDistr();
    qreal mean = 0;
    int n = in.size();
    qreal So = 1.0/(12*n);

    if (alpha == -1) {
        for (int i=0; i<n; i++)
            mean += in.at(i)/n;

        distr->setParameter(ALPHA, mean/beta);
    }
    else
        distr->setParameter(ALPHA, alpha);
    distr->setParameter(BETA, beta);

//    if (in.isEmpty())
//        in = distr->getSortedData();
//    else {
//        qSort(in);
//    }

    qSort (in);

    for (int i=1; i<=n; i++) {
        qreal delta;
//        delta = distr->cdf(in.at(i-1));
        delta = distr->AbstractDistribution::cdf(in.at(i-1)) - (2.0*i - 1.0)/(2*n);
        So += pow(delta, 2);
    }

//    qDebug () << So;

    delete distr;

    return So;
}

qreal GammaDistr::defineAlphaFromBeta(qreal beta, const array1 &in)
{
    if (in.size() < 2)
        return 0;

    qreal mean=0;
    for (int i=0; i<in.size(); i++) {
        mean += in.at(i);
    }
    mean /= in.size();

    return mean/beta;
}

qreal GammaDistr::defineBetaFromAlpha(qreal alpha, const array1 &in)
{
    if (in.size() < 2)
        return 0;

    qreal mean=0;
    for (int i=0; i<in.size(); i++) {
        mean += in.at(i);
    }
    mean /= in.size();

    return mean/alpha;
}

void GammaDistr::defineParameters(const array1 &in, array1 &pars)
{
    pars.resize(2);

    qreal mean=0, s2=0;

    for (int i=0; i<in.size(); i++) {
        mean += in.at(i);
    }

    mean /= in.size();

    for (int i=0; i<in.size(); i++) {
        s2 += pow(in.at(i) - mean, 2);
    }

    s2 /= in.size()-1;


    qreal beta0 = 0.1;
//    beta0 = getRoot(GammaEq, in, beta0);
////    qDebug () << "Eq Solved!";

//    qreal alpha0 = 0;
//    for (int i=0; i<in.size(); i++) {
//        alpha0 += in.at(i)/(in.size()*beta0);
//    }
//    qreal alpha0 = mean/beta0;

    qreal alpha0 = mean*mean/s2/*-1*/;
    beta0  = s2/mean;
    pars[BETA] = beta0;
    pars[ALPHA] = alpha0;
    pars = estimate(in, pars); return ;

    pars = findMDPars(in, pars); return;
}



qreal GammaDistr::GammaFunc(qreal z)
{
    return alglib::gammafunction(z);
}

qreal GammaDistr::gammaFunc(qreal a, qreal z)
{
    //! Неполная гамма-функция
    //! a - параметр
    //! z - переменная

    return alglib::incompletegamma(a, z)*alglib::gammafunction(a);
}

void GammaDistr::getMirroredData()
{

}

qreal GammaDistr::pdf(qreal x, const array1 &tmpPars)
{
    qreal out;
    if (x<0)
        return 0;
    out =  pow(x, tmpPars[ALPHA]-1)*exp(-x/tmpPars[BETA]);
    out /= GammaFunc(tmpPars[ALPHA])*pow(tmpPars[BETA], tmpPars[ALPHA]);
    return out;
}

//qreal GammaDistr::rnd()
//{
//    qreal y = rng.uniform();
//    return getRoot(cpd, pars, 2, y);
//}

// function b_w(theta) for Algorithm 2.1
inline qreal GammaDistr::convertl(qreal theta)
{
    qreal bw;
    qreal alpha = pars[ALPHA]/*-1*/;
    // All the numbers in the function are obtained from Table 1
    // The constant c1, c2, and c3 are used to declare the numbers for range of theta
    // The constant h1, h2, h3, h4, h5, and h6 are used to declare the numbers used to calculate b_w(theta)

    qreal c1=1.76421669;
    qreal c2=0.52122324;
    qreal c3=0.20931492;
    qreal h1=-0.04806589;
    qreal h2=-0.08476353;
    qreal h3=0.10147534;
    qreal h4=-0.13546023;
    qreal h5=0.12789964;
    qreal h6=-0.30685282;

    if(theta>=c1)
    {
        bw=h1;
        return bw;
    }
    else if(theta>=c2 && theta<c1)
    {
        bw=h2*theta+h3;
        return bw;
    }
    else if(theta>=c3 && theta<c2)
    {
        bw=h4*theta+h5;
        return bw;
    }
    else
    {
        bw=h6-theta/2+alpha/2;
        return bw;
    }
}

// function b_s(theta) for Algorithm 2.1
inline qreal GammaDistr::convertu(qreal theta)
{
    qreal bs;
    // All the numbers in the fuction are obtained from Table 2
    // The constant c1 and c2 are used to declare the numbers for range of theta
    // The constant h1,h2, h3, h4, and h5 are used to declare the numbers used to calculate b_s(theta)

    qreal c1=1.44893155;
    qreal c2=-3.33318991;
    qreal h1=-0.15342641;
    qreal h2=0.12465180;
    qreal h3=-0.33403833;
    qreal h4=0.30625300;
    qreal h5=0.27127295;

    if(theta>c1)
    {
        bs=h1;
        return bs;
    }
    else if(theta>c2 && theta<=c1)
    {
        bs=h2*theta+h3;
        return bs;
    }
    else
    {
        bs=h4*theta+h5;
        return bs;
    }
}

// Algorithm 2.1, all shape parameter pars[ALPHA]>0
qreal GammaDistr::rnd(/*int k, qreal *a*/)
{
    qreal alpha = pars[ALPHA]/*-1*/;


    qreal theta=log(alpha);
    qreal Bmin=-exp(convertl(theta));
    qreal Bmax=exp(convertu(theta));
    qreal u;
    qreal v;
    qreal c=sqrt(alpha);
    qreal t,t1;
//    int i=0;

    // Check pars[ALPHA] condition
    if(alpha>=0.01)
    {
        do
        {
            // Generate Uniform[0,1] and assign the number to variable u
            u=rng.uniform();
            // Generate Uniform[Bmin, Bmax] and assign it to v
            v=Bmin+rng.uniform()*(Bmax+(-Bmin));
            t=v/u;
            t1=exp(t/c+theta);
            // Check acceptance condition
            if(2*log(u)<=c*t-t1+alpha)
            {
                return t1*pars[BETA];
//                a[i]=t1;
//                i=i+1;
            }
        }while(true/*i<k*/);
    }
    // Generate random numbers on log scale when pars[ALPHA] is <0.01
    else if(alpha<0.01)
    {
        do
        {
            u=rng.uniform();
            v=Bmin+rng.uniform()*(Bmax+(-Bmin));
            t=v/u;
            t1=exp(t/c+theta);
            // Check acceptance condition
            if(2*log(u)<=c*t-t1+alpha)
            {
                return ((t/c)+theta)*pars[BETA];
//                a[i]=(t/c)+theta;
//                i=i+1;
            }
        }while(true/*i<k*/);
    }

    return 0;
}

bool GammaDistr::verify(const array1 &tmpPars)
{
    return tmpPars.at(BETA) > 0 && tmpPars.at(ALPHA) > 0;
}
