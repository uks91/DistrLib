#ifndef ABSTRACTDISTRIBUTION_H
#define ABSTRACTDISTRIBUTION_H

//#include <QtCore>
#include "../functions.h"
#include "../rng.h"



class AbstractDistribution
{
public:
    AbstractDistribution();
    AbstractDistribution(const array1 &in, const array1 *pars0=0);
    ~AbstractDistribution();
    enum Method {Kolmogorov, omega_2, Omega_2};
    enum Estimation {MLH, MD};

            qreal   cdf             (qreal x) {return cdf(x, pars);}

            array1  estimate        (array1 in, array1 pars0);
            array1  generateSeries  (int n);
            array1  generateSeries  ();
            array1  getData         () {return data;}
            qreal   getMinOrMax     () {return minOrMax;}
            array1  getPars         () {return pars;}
            array1  getPars         (const array1 &input);
            array1  getSortedData   () {return sortedData;}
            array1  findMDPars      (array1 in, array1 pars0);
//            array1  findMDPars      (array1 in, array1 pars0, qreal &minValue);
            array1  findMLHPars     (array1 in, array1 pars0/*, qreal &maxValue*/);
            qreal   findPValue      ();
            qreal   findPValue      (array1 in, qreal &criteria, int iters=1E4, Method method=omega_2);
            qreal   fitCriteria     (array1 in, Method method=omega_2) {return fitCriteria(in, pars, method);}
            qreal   fitCriteria     (array1 in, const array1 &tmpPars, Method method=omega_2, bool chkPars = false);
    virtual QString name            () = 0; //{return "Abstract Distribution";}
            qreal   pdf             (qreal x) {return pdf(x, pars);}
    virtual qreal   rnd             () = 0;
            void    setEmpData      (const array1 &in);
            bool    setParameter    (int id, qreal value);
            bool    setParameter    (const array1 &in);
            int     numGenerations  () {return rng.numGenerations;}

    static  Estimation est;

//    static  int count;

protected:
            qreal   cdf             (qreal x, const array1 &tmpPars);
    virtual qreal   cdfVerified     (qreal x, const array1 &tmpPars) = 0;
//            qreal   pdf             (qreal x, const array1 &tmpPars);
    virtual qreal   pdf/*Verifird*/     (qreal x, const array1 &tmpPars) = 0;
    virtual void    defineParameters(const array1 &d, array1 &p) = 0;
            void    defineParameters() {defineParameters(data, pars);}
            qreal   likelihood   (const array1 &in, const array1 &tmpPars);
    virtual bool    verify       (const array1 &tmpPars) = 0;

            array1  data;  //эмпирические значения, по которым подбираются параметры распределения
            array1  sortedData;
            array1  pars;  //параметры рпспределения
            bool    preDefined;

    static  RNG     rng;

            qreal   minOrMax;    //min значение критерия или max значение multiPDF

//private:
    qreal  explore    (array1 in, array1 h, array1 x1, array1 &x2);
    array1 HookJeevs  (array1 in, array1 pars0);
    qreal  fitCriteria(array1 in, const array1 &tmpPars, qreal x, int ind);
    void   find2DMin  (array1 in, array1 &pars0, int ind);
    array1 vectorSum  (array1 in1, array1 in2, qreal m1=1, qreal m2=1);
    array1 vectorMult (array1 in, qreal multiplicator);
    array1 randomUnitVector(array1 p/*=array1 ()*/);
};

#endif // ABSTRACTDISTRIBUTION_H

//#define MD
//#define MOMENT
