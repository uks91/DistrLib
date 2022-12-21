#include "mathematicadistr.h"

MathematicaDistr::MathematicaDistr(const QString &name, int parsCount) : AbstractDistribution()
{
    this->distrName = name;
    pars.resize(parsCount);
    preDefine = true;
    proc.start("\"c:/Program Files/Wolfram Research/Mathematica/9.0/math.exe\" -noprompt -noinit"); // -noprompt
    /*qDebug () << */proc.waitForStarted();

}

MathematicaDistr::MathematicaDistr(const QString &name, int parsCount, const array1 &in, array1 *pars0) : AbstractDistribution(in)
{
    this->distrName = name;
    pars.resize(parsCount);
    if (pars0 == 0)
        preDefine = true;
    else {
        if (pars0->size() == parsCount) {
            preDefine = false;
            pars = *pars0;
        }
        else
            preDefine = true;
    }
    proc.start("\"c:/Program Files/Wolfram Research/Mathematica/9.0/math.exe\" -noprompt -noinit"); // -noprompt
    /*qDebug () << */proc.waitForStarted();
    AbstractDistribution::defineParameters();
}

qreal MathematicaDistr::cdfVerified(qreal x, const array1 &tmpPars)
{
    return cmd2Real(QString("CDF[%1, %2]\n").arg(prepareDistrString(tmpPars)).arg(x));
}

array1 MathematicaDistr::cdfArray(array1 x/*, const array1 &tmpPars*/)
{
    return cmd2array( QString("CDF[%1, %2]\n").arg(prepareDistrString(pars)).arg(arrayToString(x)) );
}

void MathematicaDistr::defineParameters(const array1 &in, array1 &pars)
{
    return;

    if (preDefine) {
        QString data = arrayToString(in);
        QString p = parsToString();
        QString command = QString ("{%1}/.FindDistributionParameters[%2, %3[%1]]\n").arg(p).arg(data).arg(name());
//        qDebug () << command;
        pars = cmd2array(command);
//        qDebug () << "PreDefine: " << pars;
//        return;

    }
//    qDebug () << "PreDefine: " << pars;

    int N = 1000, M = 20;
    qreal step = 1, beta=0.5, e=1E-3, stepVariance=1;
//    int N = 1000, M = 20;
//    qreal step = 0.5, beta=0.5, e=1E-3;
    array1 x0, x1, ksi;
    qreal basestep = step;
    x0 = pars;
    x1 = x0;


//    if (name()=="Meixner") {
//        step = 1;
//        e = 3E-3;
//    }

    for (int i=0; i<N; i++) {
        qreal min = fitCriteriaMath(in, x0, omega_2, true);
        bool chgCenter = false;

        for (int j=0; j<M; j++) {
            ksi = vectorMult(randomUnitVector(pars/*, pars0*/), step);
            array1 tmp = vectorSum(x0, ksi);
            qreal S = fitCriteriaMath(in, tmp, omega_2, true);
            if (S < min) {
                chgCenter = true;
                min = S;
                x1 = tmp;
            }
        }

        if (chgCenter) {
            x0 = x1;
            if (stepVariance < 2)
                stepVariance += 0.05;
            step = basestep*stepVariance;
        }
        else {
            step *= beta;
            basestep = step;
            stepVariance = 1;
        }

        if (step < e) {
            pars = x0;
//            qDebug () << QString("Solved at %1 steps").arg(i+1);
            return;
        }
//        qDebug () << i << step;
    }
//    pars = x0;
}

qreal MathematicaDistr::fitCriteriaMath(array1 in, const array1 &tmpPars, Method method, bool chkPars)
{
    if (chkPars) {
        for (int i=0; i<tmpPars.size(); i++) {
            if (tmpPars.at(i) < 0)
                return 1E10;
        }
    }
    QString m = "CramerVonMisesTest";
    QString data = arrayToString(in);

    if (method == Omega_2)
        m = "AndersonDarling";
    bool b;
    qreal tmp = cmd2Real( QString("%1[%2, %3, \"TestStatistic\"]\n").arg(m).arg(data).arg(prepareDistrString(tmpPars)), &b );
    if (b)
        return tmp;
    else
        return 1E10;
}

array1 MathematicaDistr::generateSeries(int n)
{
    return cmd2array( QString("RandomVariate[%1, %2]\n").arg(prepareDistrString()).arg(n) );
}

qreal MathematicaDistr::pdf(qreal x, const array1 &tmpPars)
{
   return cmd2Real(QString("PDF[%1, %2]\n").arg(prepareDistrString(tmpPars)).arg(x));
}

qreal MathematicaDistr::rnd()
{
//    qDebug () << cmd2array( QString("RandomVariate[%1, 4]\n").arg(prepareDistrString()) );
    return cmd2Real( QString("RandomVariate[%1]\n").arg(prepareDistrString()) );
}

QString MathematicaDistr::arrayToString(const array1 &in)
{
    QString data = "{";
    for (int i=0; i<in.size(); i++) {
        data += QString().number(in.at(i));
        if (i<in.size()-1)
            data += ", ";
    }
    data += "}";
    return data;
}

array1 MathematicaDistr::cmd2array(const QString &command, bool *b)
{
    array1 out;
    QStringList list = sendCommand(command).replace("{", "").replace("}", "").split(", ");
//    qDebug () << list;
    for (int i=0; i<list.size(); i++) {
        out << list.at(i).toDouble(b);
    }
    return out;
}

qreal MathematicaDistr::cmd2Real(const QString &command, bool *b)
{
//    qDebug () << command;
    return sendCommand(command).toDouble(b);
}

QString MathematicaDistr::parsToString()
{
    QString out = "";
    for (int i=0; i<pars.size(); i++) {
        out += QString("p%1").arg(i);
        if (i!= pars.size()-1)
            out += ", ";
    }
    return out;
}

QString MathematicaDistr::parsToString(const array1 &p)
{
    QString out = "";
    for (int i=0; i<pars.size(); i++) {
        out += QString().number( p.at(i) );
        if (i!= pars.size()-1)
            out += ", ";
    }
    return out;
}

QString MathematicaDistr::prepareDistrString(array1 tmpPars)
{
    return QString("%1[%2]")
            .arg(name())
            .arg(arrayToString(tmpPars).replace("{", "").replace("}", ""));
}

QString MathematicaDistr::sendCommand(const QString &command)
{
//    proc.waitForFinished(1);
    proc.waitForReadyRead(1);
    proc.write(command.toLatin1());
//    proc.closeWriteChannel();
    proc.waitForBytesWritten();
    proc.waitForReadyRead();
    return QString(proc.readAll()).trimmed();


    return command;
}



MathematicaDistr::~MathematicaDistr()
{
    proc.close();
}
