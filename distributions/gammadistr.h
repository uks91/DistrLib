#ifndef GAMMADISTR_H
#define GAMMADISTR_H

#include "abstractdistribution.h"
//#include "alglib/specialfunctions.h"


qreal cpd (qreal x, array1 p);

class GammaDistr : public AbstractDistribution
{


public:
    //! ALPHA > 0
    enum {ALPHA, BETA, SHIFT/*, K*/};
    GammaDistr(qreal shift=0, bool isMirrored=false);
    GammaDistr(const array1 &in, const array1 *pars0=0);

    qreal cdfVerified (qreal x, const array1 &tmpPars);

    void defineAlphaFromBeta () {pars[ALPHA] = defineAlphaFromBeta(pars[BETA], data);}
    qreal defineAlphaFromBeta(qreal beta, const array1 &in);
    void defineBetaFromAlpha (){pars[BETA] = defineBetaFromAlpha(pars[ALPHA], data);}
    qreal defineBetaFromAlpha (qreal alpha, const array1 &in);
    void defineParameters(const array1 &in, array1 &pars);
    QString name () {return "GammaDistr";}
    qreal pdf (qreal x, const array1 &tmpPars); //плотность вероятности
    qreal rnd();

    qreal optimum (qreal beta, array1 in/*, qreal shift=0*/, qreal alpha=-1);

    static qreal GammaFunc (qreal z);
    static qreal gammaFunc (qreal a, qreal z);
//    static qreal gammaFunc2 (qreal a, qreal z);
//    void setEmpData(const array1 &in);
    static constexpr qreal Euler = .5772156649;


private:
    bool isMirrored;
    void getMirroredData ();
    inline qreal convertl(qreal theta);
    inline qreal convertu(qreal theta);
    bool  verify(const array1 &tmpPars);

//    void rgamma31();




    array1 mirroredData;

};

#endif // GAMMADISTR_H
