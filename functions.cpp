#include "functions.h"

void arrayToFile(const QString &fileName, const array1 &in)
{
    QFile file (fileName);
    QTextStream out (&file);
    file.open(QFile::WriteOnly);

    for (int i=0; i<in.size(); i++) {
        out << in.at(i) << endl;
    }

    file.close();
}


void arrayToFile(const QString &fileName, const array2 &in)
{
    QFile file (fileName);
    QTextStream out (&file);
    file.open(QFile::WriteOnly);

    for (int i=0; i<in.size(); i++) {
        for (int j=0; j<in.at(i).size(); j++) {
            out << in.at(i).at(j);
            if (j == in.at(i).size()-1)
                out << endl;
            else
                out << "\t";
        }
    }

    file.close();
}


array2 CDF(array1 series)
{
//    return ;
    array2 out;
    array1 tmp;

    qSort (series);
    tmp << series.first() << 0;
    out << tmp;

    for (int i=1; i<series.size(); i++) {
        tmp.clear();
        if (series.at(i) != series.at(i-1)) {
//            out << QString::number(series.at(i), 'g', 20) << "\t" << QString::number (i*1.0/series.size(), 'g', 10) << endl;
            tmp << series.at(i) << i*1.0/series.size();
            out << tmp;
        }
    }

    return out;
}

//array1 fileToArray(const QString &fileName)
//{
//    QFile file (fileName);
//    QTextStream in (&file);
//    file.open(QFile::ReadOnly);
//    array1 out;
//    while (!in.atEnd()) {
//        QString tmp;
//        in >> tmp;
//        bool b;
//        qreal n = tmp.toDouble(&b);
//        if (b)
//            out << n;
//        else
//            return out;
//    }

//    return out;
//}


array1 getMirroredData(const array1 &in, qreal margin, bool isMirrored, bool is4)
{
    array1 out;
    for (int i=0; i<in.size(); i++) {
        qreal tmp;
        if (!is4)
            tmp = in.at(i);
        else
            tmp = 90.0 + in.at(i)*1.0/10;

        if (isMirrored)
            out << margin - tmp/*in.at(i)*/;
        else
            out << tmp/*in.at(i)*/ - margin;
    }
    return out;
}

qreal getRoot(function func, array1 pars, qreal x0, qreal res)
{
    const qreal e=1e-5;

    array1 x (3);
    array1 u (2);

    x[1] = x0;  //Начальное приближение
    x[0] = x[1]*0.85;  //Второе начальное приближение

    do {

        u[0] = (*func)(x[0], pars)-res;
        u[1] = (*func)(x[1], pars)-res;

        x[2] = x[1] - u[1]*(x[1]-x[0])/(u[1]-u[0]);
        x[0] = x[1];
        x[1] = x[2];
    }
    while (fabs(/*res - */u[1])>e);

    return x[1];
}

array1 loadData(const QString &filename, int maxSize)
{
    QFile file (filename);
    if (!file.open(QFile::ReadOnly))
    {
        qDebug () << "File " << file.fileName() << "is not opened!" << file.exists();
    }
    QTextStream in (&file);
    array1 out;

    int i=0;

    while (!in.atEnd()) {
        if (i == maxSize)
            break;

        QString tmp;
        qreal number;
        bool  ok;

        in >> tmp;
        tmp.replace(",", ".");
        number = tmp.toDouble(&ok);

        if (!ok)
            break;
        out << number;

        i++;
    }
    return out;
}

QString number(qreal n, int d)
{
    int m = pow(10, d);
    n = qFloor(n*m);
    n /= 1.0*m;
    return QLocale().toString(n, 'g');
}

qreal mround (qreal x, qreal signs)
{
//    int p = qFloor(log10(x));
    if (x < pow(10, -signs))
        return x;
    qreal x2 = x*pow(10, signs);
    qreal mean = .5*(floor(x2)+ceil(x2));
    if (x2 >= mean)
        return ceil(x2)*pow(10, -signs);
    else
        return floor(x2)*pow(10, -signs);
}
