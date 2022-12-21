#ifndef MATHEMATICADISTR_H
#define MATHEMATICADISTR_H

#include "abstractdistribution.h"


class MathematicaDistr : public AbstractDistribution
{
public:
//    enum {A, B, C, D};
    MathematicaDistr(const QString &name, int parsCount);
    ~MathematicaDistr();
    MathematicaDistr(const QString &name, int parsCount, const array1 &in, array1 *pars0=0);

    array1  cdfArray        (array1 x);
    array1  generateSeries  (int n);
    QString name () {return distrName;}
    qreal   pdf             (qreal x, const array1 &tmpPars);
    qreal   rnd             ();

    qreal   fitCriteriaMath (array1 in, const array1 &tmpPars, Method method=omega_2, bool chkPars = false);

protected:
    qreal   cdfVerified     (qreal x, const array1 &tmpPars);
    void    defineParameters(const array1 &in, array1 &pars);

private:
    QString  distrName;
    QProcess proc;
    bool     preDefine;
    QString  arrayToString      (const array1 &in);
    QString  parsToString       ();
    QString  parsToString       (const array1 &p);
    QString  prepareDistrString () {return prepareDistrString(pars);}
    QString  prepareDistrString (array1 tmpPars);
    QString  sendCommand        (const QString &command);
    array1   cmd2array          (const QString &command, bool *b=0);
    qreal    cmd2Real           (const QString &command, bool *b=0);
    bool     verify(const array1 &/*tmpPars*/) {return true;}
};

#endif // MATHEMATICADISTR_H
